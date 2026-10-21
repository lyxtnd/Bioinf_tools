from sam import SamReader


def demonstrate_sam(filename: str) -> None:
    """Демонстрирует работу методов SamReader на одном SAM-файле."""
    print(f"\n{'=' * 70}")
    print(f"АНАЛИЗ SAM-ФАЙЛА: {filename}")
    print("=" * 70)

    try:
        reader = SamReader(filename)

        # 1. Заголовок и группы метаданных
        print("\n1. Заголовок файла и группы метаданных:")
        header = reader.get_header()
        print(f"Всего строк заголовка: {len(header.raw_lines)}")

        hd_group = header.get_group("HD")
        print(f"Группа @HD: {hd_group}")

        sq_group = header.get_group("SQ")
        print(f"Группа @SQ (количество последовательностей): {len(sq_group)}")
        if sq_group:
            print(f"  Пример @SQ: {sq_group[0]}")

        rg_group = header.get_group("RG")
        print(f"Группа @RG (Read Groups): {rg_group}")

        pg_group = header.get_group("PG")
        print(f"Группа @PG (Programs): {pg_group}")

        print(f"Референсный геном: {reader.get_reference_genome()}")

        # 2. Список хромосом и валидация координат
        chromosomes = reader.get_chromosomes()
        print(f"\n2. Хромосомы в файле ({len(chromosomes)} шт.): {chromosomes}")

        if chromosomes:
            test_chrom = chromosomes[0]
            print(f"Проверка координаты {test_chrom}:1000 -> {reader.validate_coordinate(test_chrom, 1000)}")
            print(f"Проверка координаты {test_chrom}:-5 -> {reader.validate_coordinate(test_chrom, -5)}")
            print(f"Проверка координаты unknown_chr:100 -> {reader.validate_coordinate('unknown_chr', 100)}")

        # 3. Подсчёт общего количества выравниваний
        total_count = reader.count_alignments()
        print(f"\n3. Общее количество выравниваний: {total_count}")

        # 4. Вывод первых выравниваний
        all_alignments = reader.read_alignments()
        if all_alignments:
            print("\nПервые 3 выравнивания:")
            for aln in all_alignments[:3]:
                print(f"  {aln}")
                print(f"    SEQ[:20]: {aln.seq[:20]} | Tags: {aln.tags}")

        # 5. Статистика 'количество выравниваний - хромосома' (pandas)
        print("\n4. Статистика выравниваний по хромосомам (pandas DataFrame):")
        df_stats = reader.count_alignments_by_chromosome()
        print(df_stats.to_string(index=False))

        # 6. Поиск в геномном интервале (аналог bedtools intersect)
        if all_alignments:
            mapped = [a for a in all_alignments if a.is_mapped]
            if mapped:
                sample_aln = mapped[0]
                chrom = sample_aln.rname
                start = max(1, sample_aln.pos - 50)
                end = sample_aln.end_pos + 50

                print(f"\n5. Поиск пересечений с интервалом {chrom}:{start}-{end} (bedtools intersect):")
                region_alignments = reader.get_alignments_in_region(chrom, start, end)
                print(f"Найдено пересечений: {len(region_alignments)}")

                for r_aln in region_alignments[:5]:
                    print(
                        f"  QNAME={r_aln.qname} | POS={r_aln.pos}-{r_aln.end_pos} | "
                        f"CIGAR={r_aln.cigar} | MAPQ={r_aln.mapq}"
                    )

    except FileNotFoundError:
        print(f"Ошибка: Файл не найден: {filename}")
    except (ValueError, OSError) as error:
        print(f"Ошибка при обработке файла: {error}")


if __name__ == "__main__":
    test_files = ["SAM/example_1.sam", "SAM/example_2.sam"]
    for sam_path in test_files:
        demonstrate_sam(sam_path)

    print("\n[OK] Демонстрация модуля SAM успешно завершена!")