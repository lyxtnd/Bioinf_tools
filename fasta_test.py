from fasta import FastaReader, get_average_length


def demonstrate_fasta(filename):
    """Демонстрирует работу методов FastaReader на одном FASTA-файле."""
    print(f"\n{'=' * 65}")
    print(f"Читаем файл: {filename}")
    print("=" * 65)

    try:
        reader = FastaReader(filename)

        #выводим первые записи
        for seq in reader.read():
            print(f"Заголовок: {seq.header}")
            print(f"Длина: {len(seq)}")
            print(f"Тип: {seq.alphabet}")
            print()

        #средняя длина — через твою функцию
        reader2 = FastaReader(filename)
        print(f"Средняя длина: {get_average_length(reader2):.2f}")

    except FileNotFoundError:
        print(f"Ошибка: Файл не найден: {filename}")
    except (ValueError, OSError) as error:
        print(f"Ошибка при обработке файла: {error}")


if __name__ == "__main__":
    test_files = ["FASTA/example_1.fasta", "FASTA/example_2.fasta"]
    for fasta_path in test_files:
        demonstrate_fasta(fasta_path)

    print("\n[OK] Демонстрация модуля FASTA успешно завершена!")