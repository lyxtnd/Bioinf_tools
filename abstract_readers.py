from abc import ABC, abstractmethod #позволяет создать абстрактный класс
from typing import Iterator


class Reader(ABC):
    """Абстрактный базовый класс для чтения биологических файлов.
    Отвечает за общую для всех форматов работу с файлом: хранит имя файла и закрывает его.
    Как именно разбирать записи, решают наследники в методе read.
    """

    def __init__(self, filename: str):
        """
        Создаёт объект для чтения файла.

        Args:
            filename: Путь к файлу, который нужно прочитать.
        """
        self.filename = filename
        self.file = None

    @abstractmethod #требует, чтобы конкретный класс реализовал этот метод
    def read(self) -> Iterator:
        """
        Читает файл и возвращает записи по одной.

        Этот метод должен быть реализован в дочернем классе.
        Iterator позволяет обрабатывать большие файлы постепенно.
        """
        pass

    def close(self) -> None:
        """Закрывает файл, если он был открыт."""
        if self.file is not None and not self.file.closed:
            self.file.close()

    @abstractmethod
    def _parse_line(self, line: str):
        """
        Обрабатывает строку файла и преобразует её в запись.

        Этот метод должен быть реализован в дочернем классе.

        Args:
            line: Строка из файла.

        Returns:
            Запись, полученная из строки.
        """
        pass


class GenomicDataReader(Reader):
    """
    Абстрактный класс для чтения геномных данных.

    Определяет методы для получения информации о хромосомах,
    референсном геноме и проверки геномных координат.
    Конкретная реализация зависит от формата файла.
    """

    @abstractmethod
    def get_chromosomes(self) -> list[str]:
        """
        Возвращает список названий хромосом, представленных в файле.

        Returns:
            list[str]: Список названий хромосом.
        """
        pass

    @abstractmethod
    def get_reference_genome(self) -> str:
        """
        Возвращает информацию о референсном геноме,
        используемом в файле, если она доступна.

        Returns:
            str: Название или идентификатор референсного генома.
        """
        pass

    @abstractmethod
    def validate_coordinate(self, chrom: str, pos: int) -> bool:
        """
        Проверяет корректность геномной координаты.

        Args:
            chrom (str): Название хромосомы.
            pos (int): Позиция на хромосоме.

        Returns:
            bool: True, если координата корректна, иначе False.
        """
        pass
