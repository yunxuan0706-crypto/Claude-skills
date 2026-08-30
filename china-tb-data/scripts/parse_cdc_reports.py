#!/usr/bin/env python3
"""
解析CDC月度法定传染病疫情报告HTML，提取肺结核数据
=================================================
前置: 先运行 download_all_sources.py 获取月报HTML文件

依赖: pip install beautifulsoup4 pandas lxml
"""

import os
import re
import pandas as pd
from pathlib import Path
from bs4 import BeautifulSoup

BASE_DIR = Path(__file__).parent.parent
REPORTS_DIR = BASE_DIR / "raw" / "channel1" / "cdc_monthly_reports"
OUTPUT_DIR = BASE_DIR / "processed"


def extract_tb_from_report(html_path):
    """从单个月报HTML中提取肺结核发病/死亡数据"""
    with open(html_path, 'r', encoding='utf-8') as f:
        soup = BeautifulSoup(f.read(), 'html.parser')

    text = soup.get_text()

    # 尝试从表格提取
    for table in soup.find_all('table'):
        for row in table.find_all('tr'):
            cells = [td.get_text(strip=True) for td in row.find_all(['td', 'th'])]
            for i, cell in enumerate(cells):
                if '肺结核' in cell or '结核' in cell:
                    numbers = []
                    for c in cells[i+1:]:
                        cleaned = re.sub(r'[^\d.]', '', c)
                        if cleaned:
                            try:
                                numbers.append(float(cleaned))
                            except ValueError:
                                pass
                    if len(numbers) >= 2:
                        return {'cases': int(numbers[0]), 'deaths': int(numbers[1])}
                    elif len(numbers) == 1:
                        return {'cases': int(numbers[0]), 'deaths': None}

    # 备选: 从正文提取
    patterns = [
        r'肺结核[^\d]*?(\d[\d,]+)\s*例[^\d]*?死亡\s*(\d[\d,]+)',
        r'肺结核[^\d]*?报告发病\s*(\d[\d,]+)',
        r'结核[^\d]*?(\d[\d,]+)\s*例',
    ]

    for pat in patterns:
        m = re.search(pat, text)
        if m:
            groups = m.groups()
            cases = int(groups[0].replace(',', ''))
            deaths = int(groups[1].replace(',', '')) if len(groups) > 1 else None
            return {'cases': cases, 'deaths': deaths}

    return None


def extract_date_from_filename(fname):
    """从文件名提取年月"""
    # Try patterns like 2024年12月, 2024_12, etc.
    m = re.search(r'(\d{4})\D*?(\d{1,2})\D*月', fname)
    if m:
        return int(m.group(1)), int(m.group(2))
    m = re.search(r'(\d{4})\D+(\d{1,2})', fname)
    if m:
        return int(m.group(1)), int(m.group(2))
    return None, None


def main():
    OUTPUT_DIR.mkdir(exist_ok=True)

    if not REPORTS_DIR.exists():
        print(f"月报目录不存在: {REPORTS_DIR}")
        print("请先运行 download_all_sources.py")
        return

    results = []
    for html_file in sorted(REPORTS_DIR.glob("*.html")):
        year, month = extract_date_from_filename(html_file.name)
        data = extract_tb_from_report(html_file)

        if data and year:
            results.append({
                'year': year,
                'month': month,
                'tb_reported_cases': data.get('cases'),
                'tb_reported_deaths': data.get('deaths'),
                'source_file': html_file.name
            })
            print(f"  {year}-{month:02d}: cases={data.get('cases')}, deaths={data.get('deaths')}")
        else:
            print(f"  [skip] {html_file.name} (无法提取)")

    if results:
        df = pd.DataFrame(results)
        df = df.sort_values(['year', 'month'])
        out_path = OUTPUT_DIR / "06_china_national_cdc_monthly_tb_reports.csv"
        df.to_csv(out_path, index=False)
        print(f"\n保存到 {out_path}: {len(df)} 条记录")
    else:
        print("\n未提取到数据，请检查HTML文件内容")


if __name__ == "__main__":
    main()
