import re
import pandas as pd
from abstract_readers import GenomicDataReader


class Alignment:
    """Представляет отдельное выравнивание прочтения из SAM-файла.

    Хранит информацию об имени прочтения, флаге, хромосоме, позиции,
    CIGAR-строке, качестве и дополнительных опциональных тегах.

    Attributes:
        qname (str): Имя прочтения (Query template NAME).
        flag (int): Битовый флаг выравнивания.
        rname (str): Название хромосомы/референса.
        pos (int): Геномная позиция начала выравнивания (1-based).
        mapq (int): Качество выравнивания (Mapping Quality).
        cigar (str): CIGAR-строка операций выравнивания.
        rnext (str): Хромосома следующего прочтения в паре.
        pnext (int): Позиция следующего прочтения в паре.
        tlen (int): Наблюдаемая длина фрагмента (Template Length).
        seq (str): Нуклеотидная последовательность прочтения.
        qual (str): Строка качества Phred.
        tags (dict[str, str]): Дополнительные теги выравнивания.
    """

    def __init__(
        self,
        qname: str,
        flag: int,
        rname: str,
        pos: int,
        mapq: int,
        cigar: str,
        rnext: str,
        pnext: int,
        tlen: int,
        seq: str,
        qual: str,
        *optional_tags: str
    ):
        """Инициализирует объект выравнивания.

        Args:
            qname: Имя прочтения.
            flag: Битовый флаг.
            rname: Название хромосомы.
            pos: Позиция начала выравнивания.
            mapq: Оценка качества картирования.
            cigar: CIGAR-строка.
            rnext: Хромосома пары.
            pnext: Позиция пары.
            tlen: Длина вставки.
            seq: Последовательность прочтения.
            qual: Качество нуклеотидов.
            *optional_tags: Необязательные поля тегов формата TAG:TYPE:VALUE.
        """
        self.qname = qname
        self.flag = flag
        self.rname = rname
        self.pos = pos
        self.mapq = mapq
        self.cigar = cigar
        self.rnext = rnext
        self.pnext = pnext
        self.tlen = tlen
        self.seq = seq
        self.qual = qual

        self.tags: dict[str, str] = {}
        for tag_entry in optional_tags:
            parts = tag_entry.split(":", 2)
            if len(parts) == 3:
                self.tags[parts[0]] = parts[2]

    @property
    def is_mapped(self) -> bool:
        """Проверяет, картировано ли прочтение (флаг 0x4)."""
        return (self.flag & 4) == 0 and self.rname != "*"

    @property
    def end_pos(self) -> int:
        """Вычисляет конечную координату выравнивания на референсе по CIGAR."""
        if not self.is_mapped or self.pos <= 0:
            return self.pos

        if self.cigar == "*" or not self.cigar:
            return self.pos + len(self.seq) - 1

        ref_length = 0
        operations = re.findall(r"(\d+)([MIDNSHP=X])", self.cigar)
        for length, op in operations:
            if op in ("M", "D", "N", "=", "X"):
                ref_length += int(length)

        return self.pos + max(0, ref_length - 1)

    def __str__(self) -> str:
        return (
            f"Alignment(QNAME='{self.qname}', RNAME='{self.rname}', "
            f"POS={self.pos}-{self.end_pos}, MAPQ={self.mapq}, CIGAR='{self.cigar}')"
        )


class SamHeader:
    """Хранит информацию из заголовка SAM-файла.

    Attributes:
        raw_lines (list[str]): Исходные строки заголовка, начинающиеся с '@'.
        headers (dict[str, list[dict[str, str]]]): Разобранные группы метаданных
            по тегам (HD, SQ, RG, PG, CO).
        reference (str | None): Информация о референсном геноме или None.
    """

    def __init__(self, raw_lines: list[str]):
        """Инициализирует объект заголовка SAM-файла.

        Args:
            raw_lines: Список исходных строк заголовка.
        """
        self.raw_lines = raw_lines
        self.headers: dict[str, list[dict[str, str]]] = {
            "HD": [],
            "SQ": [],
            "RG": [],
            "PG": [],
            "CO": []
        }
        self.reference: str | None = None
        self._parse_headers()

    def _parse_headers(self) -> None:
        """Разбирает строки заголовка на группы по ключевым тегам."""
        for line in self.raw_lines:
            line_str = line.strip()
            if not line_str.startswith("@"):
                continue

            parts = line_str.split("\t")
            tag = parts[0][1:]

            if tag == "CO":
                comment_text = line_str[4:] if len(line_str) > 4 else ""
                self.headers.setdefault("CO", []).append({"text": comment_text})
                if not self.reference and ("reference" in comment_text.lower() or "assembly" in comment_text.lower()):
                    self.reference = comment_text
                continue

            field_dict: dict[str, str] = {}
            for item in parts[1:]:
                if ":" in item:
                    k, v = item.split(":", 1)
                    field_dict[k] = v

            self.headers.setdefault(tag, []).append(field_dict)

            if tag == "SQ" and not self.reference and "AS" in field_dict:
                self.reference = field_dict["AS"]

    def get_group(self, tag: str) -> list[dict[str, str]]:
        """Возвращает информацию по заданной группе заголовков.

        Args:
            tag: Название группы, например 'PG', '@PG', 'RG', '@RG'.

        Returns:
            Список словарей с атрибутами для данной группы заголовков.
        """
        clean_tag = tag.lstrip("@").upper()
        return self.headers.get(clean_tag, [])


class SamReader(GenomicDataReader):
    """Читает SAM-файл и предоставляет методы для работы с выравниваниями.

    Позволяет разбирать отдельные строки выравниваний, читать записи
    последовательно через генератор, получать заголовок с метаданными,
    считать статистику по хромосомам и искать выравнивания в регионах.
    """

    def _parse_line(self, line: str) -> Alignment | None:
        """Разбирает строку SAM и создаёт объект Alignment.

        Args:
            line: Строка файла, содержащая информацию о выравнивании.

        Returns:
            Объект Alignment при успешном разборе или None,
            если строка содержит недостаточно полей либо некорректные типы.
        """
        parts = line.rstrip("\r\n").split("\t")
        if len(parts) < 11:
            return None

        try:
            flag = int(parts[1])
            pos = int(parts[3])
            mapq = int(parts[4])
            pnext = int(parts[7])
            tlen = int(parts[8])
        except ValueError:
            return None

        return Alignment(
            parts[0],
            flag,
            parts[2],
            pos,
            mapq,
            parts[5],
            parts[6],
            pnext,
            tlen,
            parts[9],
            parts[10],
            *parts[11:]
        )

    def read(self):
        """Последовательно читает выравнивания из SAM-файла.

        Пропускает пустые строки и строки заголовка.

        Yields:
            Alignment: Очередное разобранное выравнивание.
        """
        with open(self.filename, "r", encoding="utf-8", errors="replace") as file:
            for line in file:
                if not line.strip() or line.startswith("@"):
                    continue

                alignment = self._