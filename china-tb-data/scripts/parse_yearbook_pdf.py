#!/usr/bin/env python3
"""
从中国卫生健康统计年鉴PDF提取各省肺结核发病率表格
================================================
依赖: pip install tabula-py pandas
前置: 需要安装 Java Runtime Environment (JRE)
"""

import pandas as pd
from pathlib import Path

try:
    import tabula
except ImportError:
    print("需要安装: pip install tabula-py")
    print("还需要安装 Java: https://adoptium.net/")
    exit(1)

BASE_DIR = Path(__file__).parent.parent
PDF_DIR = BASE_DIR / "raw" / "channel1" / "nhc_yearbook"
OUTPUT_DIR = BASE_DIR / "processed"


def extract_tb_tables(pdf_path):
    """
    提取年鉴中甲乙类法定报告传���病发病率表
    年鉴中通常在 '卫生防疫' 章节
    """
    print(f"正在扫描PDF: {pdf_path}")
    print("这可能需要几分钟...")

    # 尝试提取所有表格
    all_tables = tabula.read_pdf(
        str(pdf_path),
        pages='all',
        multiple_tables=True,
        pandas_options={'header': None}
    )

    print(f"共找到 {len(all_tables)} 个表格")

    tb_tables = []
    for i, table in enumerate(all_tables):
        text = table.to_string()
        if '肺结核' in text or '结核' in text:
            tb_tables.append((i, table))
            print(f"  表格 {i}: 包含肺结核数据 ({table.shape})")

    if not tb_tables:
        print("\n未在自动提取的表格中找到肺结核数据")
        print("建议:")
        print("1. 用 Adobe Reader 打开PDF，找到相关表格页码")
        print("2. 使用 tabula.read_pdf(pdf, pages='100-110') 指定页码范围")
        print("3. 或使用 Tabula GUI工具: https://tabula.technology/")
        return None

    return tb_tables


def main():
    OUTPUT_DIR.mkdir(exist_ok=True)

    pdf_files = list(PDF_DIR.glob("*.pdf"))
    if not pdf_files:
        print(f"未找到PDF文件: {PDF_DIR}")
        print("请先运行 download_all_sources.py")
        return

    for pdf_path in pdf_files:
        tables = extract_tb_tables(pdf_path)
        if tables:
            for idx, (table_num, df) in enumerate(tables):
                out_path = OUTPUT_DIR / f"yearbook_{pdf_path.stem}_table{table_num}.csv"
                df.to_csv(out_path, index=False)
                print(f"  保存: {out_path}")


if __name__ == "__main__":
    main()
