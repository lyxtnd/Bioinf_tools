from typing import Iterator

from abstract_readers import Reader


class Seq:
    """
    Класс для работы с биологической последовательностью
    Хранит саму последовательность и информацию о ней (заголовок)
    Умеет красиво печататься в формате FASTA, считать длину и определять, белок это или нуклеотид
    """   
    def __init__(self, header, sequence):
        """
        Создаем новую последовательность
        
        header — это заголовок (мета-информация, название)
        sequence — это сама последовательность (буквы)
        """
        self.__header = header
        self.__sequence = sequence

#геттеры
    @property
    def header(self):
        """Возвращает заголовок последовательности"""
        return self.__header

    @property
    def sequence(self):
        """Возвращает саму последовательность"""
        return self.__sequence
#магические методы
    def __len__(self):
        """Возвращает длину последовательности (чтобы работала функция len())"""
        return len(self.__sequence)

    def __str__(self):
        """
        Красиво печатает последовательность в формате FASTA
        То есть сначала заголовок со знаком >, потом сама последовательность
        """
        return f">{self.__header}\n{self.__sequence}"

    def __repr__(self):
        """
        Техническое представление объекта (для отладки)
        Показывает заголовок и длину
        """
        return f"Seq(header='{self.__header}', length={len(self)})"
#свойство для определения типа последовательности
    @property 
    def alphabet(self):
        """
        Определяет, какой алфавит у последовательности: белковый или нуклеотидный
        Работает по стандарту IUPAC
        Возвращает 'protein', 'nucleotide' или 'unknown'
        """
        #множество букв, которые есть только в белках
        protein_only = set("EFILPQZJ")
        #множество всех букв, которые есть в нуклеотидах IUPAC
        nucleotide_IUPAC = set("ACGTURYSWKMBDHVN")
        #превращаем последовательность в множество уникальных букв
        seq_letters = set(self.__sequence)
        #если пересечение не пустое, есть буквы белков - это белок
        if seq_letters & protein_only:
            return "protein"
        #если все буквы есть входят в набор нуклеотидов IUPAC - это нуклеотид
        elif seq_letters <= nucleotide_IUPAC:
            return "nucleotide"
        #если ничего не подошло
        else:
            return "unknown"


class FastaReader(Reader):
    """Класс для чтения файлов в формате FASTA."""
    def __init__(self, filepath):
        #вызываем конструктор родителя (Reader), чтобы сохранить filepath
        super().__init__(filepath) 

    def is_valid_fasta(self):
            """
            Проверяет, что файл действительно в формате FASTA
            Для этого смотрит, начинается ли первая непустая строка со знака >
            Возвращает True, если формат правильный, и False, если нет
            """
            with open(self.filename, 'r', encoding='utf-8') as file:
                for line in file:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith('>'):
                        return True
            return False 

    def _parse_line(self, line):
            """
            Обрабатывает одну строку файла и преобразует её в запись
            Так как одна запись FASTA занимает несколько строк (заголовок + последовательность),
            этот метод возвращает не готовую запись, а тип строки и её содержимое
            
            Возвращает кортеж ('header', текст) или ('seq', кусок), либо None для пустых строк
            """
            line = line.strip()
            if not line:
                return None
            if line.startswith('>'):
                return ('header', line[1:])
            return ('seq', line)

    def read(self):
            """
            Генератор, который читает файл по записям
            Находит заголовок (строка со знаком >), собирает последовательность и отдает готовый объект Seq
            
            Используется yield, чтобы не грузить весь файл в память сразу — это важно для больших файлов
            
            Если файл не в формате FASTA — выдает ошибку ValueError
            """
            #проверяем, что файл в нужном формате
            if not self.is_valid_fasta():
                raise ValueError("Файл не соответствует формату FASTA")
    
            #создаем переменные для хранения текущей записи
            header = None    #пока заголовка нет
            seq_chunks = []  #пока последовательность пустая
    
            with open(self.filename, 'r', encoding='utf-8') as file:
                for line in file:
                    parsed = self._parse_line(line)
                    if parsed is None:
                        continue
                    
                    kind, value = parsed
                    if kind == 'header':
                        if header is not None:
                            yield Seq(header, "".join(seq_chunks))
                        header = value
                        seq_chunks = []
                    else:
                        seq_chunks.append(value)
    
            if header is not None:
                yield Seq(header, "".join(seq_chunks))



#функция для подсчета средней длины
def get_average_length(reader):
    """
    Вычисляет среднюю длину последовательностей за один проход по генератору
    
    reader — объект класса-наследника Reader (например, FastaReader)
    Возвращает среднюю длину как float, или 0.0 если файл пустой
    """
    total_length = 0 
    count = 0

    for seq in reader.read():
        total_length += len(seq)
        count += 1
        
    if count == 0:
        return 0.0

    return total_length /  count 