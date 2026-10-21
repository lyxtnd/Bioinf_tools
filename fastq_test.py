from FASTQ import FastqRecord,FastqParser, plot_all

path=input("Введите путь к файлу:")
compressed = input("Файл сжат gzip? (y/n): ").lower() == "y"
p = FastqParser(path,compressed)
print("Записей:", p.count_sequences())
print("Средняя длина:", p.average_length())
def parser_factory():
    return FastqParser(path, compressed)
plot_all(parser_factory, "report.png")
print("Готово: report.png")