import pandas as pd
from abstract_readers import GenomicDataReader

class Variant:
    """Представляет генетический вариант из VCF-файла.

       Хранит информацию о хромосоме, позиции, аллелях, качестве,
       статусе фильтрации и дополнительных данных варианта.

       Attributes:
           chrom (str): Название хромосомы.
           pos (int): Геномная позиция варианта.
           id (str | None): Идентификатор варианта или None, если он отсутствует.
           ref (str): Референсный аллель.
           alt (str): Альтернативный аллель.
           qual (float | None): Оценка качества варианта или None.
           filter (str): Статус фильтрации варианта.
           info (list[str]): Дополнительная информация из поля INFO.
           format (str | None): Формат данных образцов из поля FORMAT.
           samples (list[str]): Данные отдельных образцов.
       """

    def __init__(self, chrom: str, pos:int, variant_id:str|None, ref:str, alt:str, qual:float|None, filter: str, info:list[str], *args):
        """Инициализирует объект генетического варианта.

               Args:
                   chrom: Название хромосомы.
                   pos: Геномная позиция варианта.
                   variant_id: Идентификатор варианта или None.
                   ref: Референсный аллель.
                   alt: Альтернативный аллель.
                   qual: Оценка качества варианта или None.
                   filter: Статус фильтрации.
                   info: Список дополнительных сведений о варианте.
                   *args: Необязательные данные, где первый элемент —
                       поле FORMAT, а остальные — данные образцов.
               """
        self.chrom = chrom
        self.pos = pos
        self.id = variant_id
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = filter
        self.info = info
        self.format = None
        self.samples = []

        if args:
            self.format = args[0]
            self.samples = list(args[1:])


class VcfHeader:
    """Хранит информацию из заголовка VCF-файла.

        Attributes:
            metadata (list[str]): Строки метаданных, начинающиеся с '##'.
            columns (list[str]): Названия столбцов заголовка.
            samples (list[str]): Названия образцов.
            file_format (str | None): Версия формата VCF или None.
            reference (str | None): Информация о референсном геноме или None.
        """
    def __init__(self, metadata: list[str], columns:list[str], samples:list[str], file_format:str|None, reference:str|None):
        """Инициализирует объект заголовка VCF-файла.

               Args:
                   metadata: Список строк метаданных.
                   columns: Список названий столбцов.
                   samples: Список названий образцов.
                   file_format: Версия формата VCF или None.
                   reference: Информация о референсном геноме или None.
               """
        self.metadata = metadata
        self.columns = columns
        self.samples = samples
        self.file_format = file_format
        self.reference = reference

class Genotype:
    """Представляет генотип из поля GT VCF-записи.

        Attributes:
            gt (str): Исходная строка генотипа.
            separator (str): Разделитель аллелей: '/' или '|'.
            gt_ls (list[int | None]): Список индексов аллелей.
                None обозначает отсутствующий аллель.
        """

    def __init__(self, gt: str):
        """Инициализирует объект генотипа и разбирает его аллели.

                Args:
                    gt: Строка генотипа, например '0/1', '1|1' или './.'.

                Raises:
                    ValueError: Если значение аллеля не является целым числом
                        и не обозначено точкой.
                """
        self.gt = gt
        self.separator = "|" if "|" in gt else "/"

        try:
            self.gt_ls = [
                int(allele) if allele != "." else None
                for allele in gt.replace("|", "/").split("/")    ]
        except ValueError:
            raise ValueError(f"Некорректный генотип: {gt}")

class VcfReader(GenomicDataReader):
    """Читает VCF-файл и предоставляет методы для работы с вариантами.
        Позволяет разбирать отдельные строки с вариантами, читать записи
        последовательно, получать список вариантов и извлекать заголовок
        файла с метаданными и информацией об образцах.
        """

    def _parse_line(self, line: str) -> Variant | None:
        """Разбирает строку VCF и создаёт объект Variant.

        Args:
            line: Строка файла, содержащая информацию о варианте.

        Returns:
            Объект Variant при успешном разборе или None,
            если строка содержит недостаточно полей либо
            некорректную позицию.
        """
        ls = line.rstrip("\n").split("\t")

        if len(ls) < 8:
            return None

        try:
            pos = int(ls[1])
        except ValueError:
            return None

        variant_id = None if ls[2] == "." else ls[2]

        try:
            qual = None if ls[5] == "." else float(ls[5])
        except ValueError:
            qual = None

        info = [] if ls[7] == "." else ls[7].split(";")

        return Variant(
            ls[0],
            pos,
            variant_id,
            ls[3],
            ls[4],
            qual,
            ls[6],
            info,
            *ls[8:]
        )

    def read(self):
        """Последовательно читает варианты из VCF-файла.
        Пропускает пустые строки и строки заголовка.

        Yields:
            Variant: Очередной разобранный генетический вариант.
        """
        with open(self.filename, "r") as file:
            for line in file:
                if not line.strip() or line.startswith("#"):
                    continue

                variant = self._parse_line(line)
                if variant is not None:
                    yield variant

    def read_variants(self) -> list[Variant]:
        """Считывает все варианты из VCF-файла в список.

        Returns:
            Список объектов Variant.
        """
        return list(self.read())

    def get_header(self) -> VcfHeader:
        """Извлекает заголовок VCF-файла.

        Собирает строки метаданных, названия столбцов,
        названия образцов, версию формата и сведения
        о референсном геноме.

        Returns:
            Объект VcfHeader с информацией из заголовка.
        """
        metadata: list[str] = []
        columns: list[str] = []

        with open(self.filename, "r") as file:
            for line in file:
                if line.startswith("##"):
                    metadata.append(line.strip())
                elif line.startswith("#CHROM"):
                    columns = line.strip().split("\t")
                elif not line.startswith("#"):
                    break

        file_format: str | None = None
        reference: str | None = None

        for line in metadata:
            if line.startswith("##fileformat="):
                file_format = line.split("=", 1)[1].strip()
            elif line.startswith("##reference="):
                reference = line.split("=", 1)[1].strip()

        samples = columns[9:] if len(columns) > 9 else []

        return VcfHeader(
            metadata,
            columns,
            samples,
            file_format,
            reference
        )

    def get_chromosomes(self) -> list[str]:
        """
        Возвращает список уникальных хромосом, в которых обнаружены варианты.

        Хромосомы собираются непосредственно из записей вариантов VCF.
        Порядок хромосом соответствует их первому появлению в файле.

        Returns:
            list[str]: Список уникальных названий хромосом,
            в которых присутствует хотя бы один вариант.
        """
        chromosomes_list = []
        chromosomes_set = set()

        for var in self.read():
            if var.chrom not in chromosomes_set:
                chromosomes_set.add(var.chrom)
                chromosomes_list.append(var.chrom)

        return chromosomes_list


    def get_reference_genome(self) -> str:
        """Возвращает информацию о референсном геноме.

        Извлекает название или описание референсного генома
        из метаданных заголовка VCF-файла.

        Returns:
            Информация о референсном геноме или сообщение,
            если она не указана в заголовке.
        """
        reference = self.get_header().reference
        if reference:
            return reference
        else:
            return 'Информация о референсном геноме не указана'

    def validate_coordinate(self, chrom: str, pos: int) -> bool:
        """
        Проверяет корректность геномной координаты.

        Проверяет, существует ли указанная хромосома среди хромосом,
        содержащих варианты, является ли позиция положительной и не
        превышает ли она длину хромосомы, если длина указана в метаданных.

        Args:
            chrom: Название хромосомы.
            pos: Геномная позиция для проверки.

        Returns:
            True, если координата считается допустимой, иначе False.
        """
        if pos < 1:
            return False

        chromosomes_list = self.get_chromosomes()
        if chrom not in chromosomes_list:
            return False

        metadata = self.get_header().metadata

        for line in metadata:
            if not line.startswith("##contig="):
                continue

            attributes = line.lstrip("##contig=<").rstrip(">").split(",")

            contig_id = None
            length_chrom = None

            for element in attributes:
                key,value = element.split("=")

                if key == "ID":
                    contig_id = value
                elif key == "length":
                    try:
                        length_chrom = int(value)
                    except ValueError:
                        return False

            if contig_id == chrom and length_chrom is not None:
                if pos > length_chrom:
                    return False
                return True

    def count_variants(self) -> int:
        """Подсчитывает общее количество вариантов в VCF-файле.

        Считывает варианты последовательно, не сохраняя их все
        одновременно в памяти.

        Returns:
            Общее количество записей с вариантами.
        """
        cnt = 0
        for _ in self.read():
            cnt += 1
        return cnt


    def filter_by_quality(self, min_qual: float) -> list[Variant]:
        """Возвращает варианты с качеством не ниже заданного порога.

        Варианты без указанной оценки качества пропускаются.

        Args:
            min_qual: Минимально допустимая оценка качества.

        Returns:
            Список вариантов, удовлетворяющих условию качества.
        """
        variant_list = []
        for var in self.read():
            if var.qual is not None:
                if var.qual >= min_qual:
                    variant_list.append(var)

        return variant_list


    def get_genotype(self, sample: str, variant: Variant) -> Genotype | None:
        """Возвращает генотип указанного образца для заданного варианта.

        Находит образец в заголовке и поле GT в формате FORMAT,
        после чего извлекает соответствующее значение из данных образца.

        Args:
            sample: Название образца.
            variant: Объект генетического варианта.

        Returns:
            Объект Genotype, если генотип найден, иначе None.
        """
        samples = self.get_header().samples

        if not samples or not variant.format:
            return None

        format_ls = variant.format.split(':')

        if 'GT' not in format_ls or sample not in samples:
            return None

        list_index = samples.index(sample)
        format_index = format_ls.index('GT')

        if list_index >= len(variant.samples):
            return None

        sample_fields = variant.samples[list_index].split(':')

        if format_index >= len(sample_fields):
            return None

        return Genotype(sample_fields[format_index])


    def count_variants_by_chromosome(self) -> pd.DataFrame:
        """Подсчитывает количество вариантов на каждой хромосоме.

        Собирает названия хромосом в таблицу pandas, группирует записи
        по хромосоме и подсчитывает количество вариантов в каждой группе.

        Returns:
            DataFrame со столбцами chrom и count, содержащими
            названия хромосом и количество вариантов соответственно.
        """
        data = []
        for var in self.read():
            data.append({'chrom': var.chrom})

        df = pd.DataFrame(data)
        result = df.groupby('chrom').size().reset_index(name='count')
        return result


    def get_variants_in_region(
        self, chrom: str, start: int, end: int
    ) -> list[Variant]:
        """Возвращает варианты, расположенные в заданном геномном интервале.

        Учитывает только варианты указанной хромосомы, позиция которых
        находится между началом и концом интервала включительно.

        Args:
            chrom: Название хромосомы.
            start: Начальная позиция интервала, начиная с 1.
            end: Конечная позиция интервала.

        Returns:
            Список вариантов, попадающих в заданный интервал.

        Raises:
            ValueError: Если начальная позиция меньше 1 или конец
                интервала расположен раньше его начала.
        """
        if start > end or start < 1:
            raise ValueError('Неверно заданы координаты')

        variants_list = []
        for var in self.read():
            if var.chrom == chrom and start <= var.pos <= end:
                variants_list.append(var)

        return variants_list



