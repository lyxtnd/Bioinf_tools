from vcf import VcfReader
def demonstrate_vcf(filename: str) -> None:
    """Демонстрирует работу методов VcfReader на одном VCF-файле."""
    print(f"\n{'=' * 65}")
    print(f"АНАЛИЗ ФАЙЛА: {filename}")
    print("=" * 65)

    try:
        reader = VcfReader(filename)

        # 1. Чтение вариантов
        variants = reader.read_variants()
        print(f"\n1. Всего вариантов: {reader.count_variants()}")

        if variants:
            print("\nПервые 5 вариантов:")
            for variant in variants[:5]:
                print(
                    f"  {variant.chrom}:{variant.pos} "
                    f"{variant.ref}>{variant.alt}, "
                    f"QUAL={variant.qual}, "
                    f"FILTER={variant.filter}"
                )
        else:
            print("В файле нет вариантов.")

        # 2. Чтение заголовка
        header = reader.get_header()
        print("\n2. Заголовок файла")
        print(f"Формат VCF: {header.file_format}")
        print(f"Референсный геном: {reader.get_reference_genome()}")
        print(f"Первые 20 колонок: {header.columns[:20]}")
        print(f"Первые 20 oбразцов: {header.samples[:20]}")
        print(f"Количество строк метаданных: {len(header.metadata)}")

        # 3. Получение хромосом
        chromosomes = reader.get_chromosomes()
        print("\n3. Хромосомы, в которых есть варианты:")
        print(chromosomes)

        # 4. Проверка координат
        if variants:
            variant = variants[0]
            print("\n4. Проверка координат")
            print(
                f"Координаты {variant.chrom}:{variant.pos}: "
                f"{reader.validate_coordinate(variant.chrom, variant.pos)}"
            )
            print(
                f"Координаты chr_UNKNOWN:0: "
                f"{reader.validate_coordinate('chr_UNKNOWN', 0)}"
            )

        # 5. Фильтрация по качеству
        print("\n5. Фильтрация вариантов по качеству (QUAL >= 30)")
        high_quality = reader.filter_by_quality(30.0)
        print(f"Найдено вариантов: {len(high_quality)}")

        for variant in high_quality[:5]:
            print(f"  {variant.chrom}:{variant.pos}, QUAL={variant.qual}")

        # 6. Получение генотипа
        print("\n6. Получение генотипа")
        if header.samples and variants:
            sample = header.samples[0]
            variant = variants[0]
            genotype = reader.get_genotype(sample, variant)

            if genotype is not None:
                print(f"Образец: {sample}")
                print(f"Вариант: {variant.chrom}:{variant.pos}")
                print(f"Генотип: {genotype.gt}")
                print(f"Разделитель аллелей: {genotype.separator}")
                print(f"Аллели: {genotype.gt_ls}")
            else:
                print("Не удалось получить генотип для выбранного варианта.")
        else:
            print("В файле нет образцов или вариантов для проверки генотипа.")

        # 7. Подсчёт вариантов по хромосомам
        print("\n7. Количество вариантов по хромосомам")
        chromosome_stats = reader.count_variants_by_chromosome()
        print(chromosome_stats.to_string(index=False))

        # 8. Поиск вариантов в регионе
        print("\n8. Поиск вариантов в регионе")
        if variants:
            first_variant = variants[0]
            chrom = first_variant.chrom
            start = max(1, first_variant.pos - 1000)
            end = first_variant.pos + 1000

            region_variants = reader.get_variants_in_region(
                chrom, start, end
            )

            print(f"Регион: {chrom}:{start}-{end}")
            print(f"Количество вариантов в регионе: {len(region_variants)}")

            for variant in region_variants[:5]:
                print(f"  {variant.chrom}:{variant.pos} {variant.ref}>{variant.alt}")
        else:
            print("Поиск региона пропущен: в файле нет вариантов.")

    except FileNotFoundError:
        print(f"Файл не найден: {filename}")
    except (ValueError, OSError) as error:
        print(f"Ошибка при обработке файла: {error}")


if __name__ == "__main__":
    files = ["VCF/example_1.vcf", "VCF/example_2.vcf", "VCF/example_3.vcf", "VCF/example_4.vcf"]
    for filename in files:
        demonstrate_vcf(filename)

    print("\nДемонстрация завершена!")

