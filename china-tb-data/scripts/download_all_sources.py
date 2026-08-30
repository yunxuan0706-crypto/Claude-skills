#!/usr/bin/env python3
"""
中国肺结核发病率数据采集脚本
=============================
在本地机器上运行此脚本以下载所有公开可获取的TB数据源。
需要网络访问中国政府网站、学术出版商和国际组织网站。

依赖: pip install requests beautifulsoup4 pandas openpyxl lxml
用法: python download_all_sources.py
"""

import os
import re
import sys
import time
import json
import requests
from pathlib import Path
from datetime import datetime

BASE_DIR = Path(__file__).parent.parent
RAW_DIR = BASE_DIR / "raw"
CHANNEL1_DIR = RAW_DIR / "channel1"
CHANNEL2_DIR = RAW_DIR / "channel2"

HEADERS = {
    "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                  "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
}

def ensure_dirs():
    for d in [CHANNEL1_DIR, CHANNEL2_DIR,
              CHANNEL1_DIR / "cdc_monthly_reports",
              CHANNEL1_DIR / "nhc_yearbook",
              CHANNEL2_DIR / "paper_supplements",
              BASE_DIR / "processed"]:
        d.mkdir(parents=True, exist_ok=True)


def download_file(url, path, desc=""):
    if path.exists():
        print(f"  [skip] {desc or path.name} (already exists)")
        return True
    print(f"  [download] {desc or path.name}")
    try:
        r = requests.get(url, headers=HEADERS, timeout=60, stream=True)
        r.raise_for_status()
        with open(path, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
        print(f"    -> saved ({path.stat().st_size:,} bytes)")
        return True
    except Exception as e:
        print(f"    -> FAILED: {e}")
        return False


# ============================================================
# 渠道一: 公共平台数据
# ============================================================

def download_cdc_monthly_reports():
    """
    采集国家疾控局/CDC月度法定传染病疫情概况
    数据来源: https://www.ndcpa.gov.cn 和 https://www.chinacdc.cn
    """
    print("\n" + "="*60)
    print("1. 国家疾控局/CDC月度传染病疫情报告")
    print("="*60)

    out_dir = CHANNEL1_DIR / "cdc_monthly_reports"

    # 国家疾控局疫情信息列表页
    print("\n从国家疾控局采集月报链接...")
    base_url = "https://www.ndcpa.gov.cn/jbkzzx/c100016/second/list.html"

    try:
        from bs4 import BeautifulSoup
        r = requests.get(base_url, headers=HEADERS, timeout=30)
        r.encoding = 'utf-8'
        soup = BeautifulSoup(r.text, 'html.parser')

        links = []
        for a in soup.find_all('a', href=True):
            text = a.get_text(strip=True)
            if '法定传染病' in text and '疫情概况' in text:
                href = a['href']
                if not href.startswith('http'):
                    href = "https://www.ndcpa.gov.cn" + href
                links.append((text, href))

        print(f"  找到 {len(links)} 个月报链接")

        for title, url in links:
            fname = re.sub(r'[^\w一-鿿]', '_', title) + ".html"
            fpath = out_dir / fname
            if not fpath.exists():
                try:
                    r2 = requests.get(url, headers=HEADERS, timeout=30)
                    r2.encoding = 'utf-8'
                    with open(fpath, 'w', encoding='utf-8') as f:
                        f.write(r2.text)
                    print(f"  [saved] {title}")
                    time.sleep(1)
                except Exception as e:
                    print(f"  [failed] {title}: {e}")
            else:
                print(f"  [skip] {title}")

    except ImportError:
        print("  需要安装 beautifulsoup4: pip install beautifulsoup4")
    except Exception as e:
        print(f"  采集失败: {e}")

    # 备选: 中国CDC结核病中心月报
    print("\n从中国CDC结核病控制中心采集月报...")
    tb_url = "https://tb.chinacdc.cn/yqjk/myyqbg/"
    try:
        from bs4 import BeautifulSoup
        r = requests.get(tb_url, headers=HEADERS, timeout=30)
        r.encoding = 'utf-8'
        soup = BeautifulSoup(r.text, 'html.parser')

        for a in soup.find_all('a', href=True):
            text = a.get_text(strip=True)
            if text and ('月' in text or '年' in text):
                href = a['href']
                if not href.startswith('http'):
                    href = "https://tb.chinacdc.cn" + href
                fname = re.sub(r'[^\w一-鿿]', '_', text) + ".html"
                fpath = out_dir / "tb_center" / fname
                fpath.parent.mkdir(exist_ok=True)
                if not fpath.exists():
                    try:
                        r2 = requests.get(href, headers=HEADERS, timeout=30)
                        r2.encoding = 'utf-8'
                        with open(fpath, 'w', encoding='utf-8') as f:
                            f.write(r2.text)
                        time.sleep(0.5)
                    except:
                        pass
    except Exception as e:
        print(f"  {e}")


def download_nhc_yearbook():
    """
    下载中国卫生健康统计年鉴PDF中的传染病数据表
    """
    print("\n" + "="*60)
    print("2. 中国卫生健康统计年鉴")
    print("="*60)

    out_dir = CHANNEL1_DIR / "nhc_yearbook"

    # 2023年鉴PDF (国家卫健委直接链接)
    yearbook_url = ("https://www.nhc.gov.cn/mohwsbwstjxxzx/tjtjnj/202501/"
                    "8193a8edda0f49df80eb5a8ef5e2547c/files/1740022743894_10341.pdf")
    download_file(yearbook_url, out_dir / "china_health_yearbook_2023.pdf",
                  "2023年中国卫生健康统计年鉴")

    print("\n  提示: 年鉴PDF中 '甲乙类法定报告传染病发病率' 表包含31省肺结核年度发病率")
    print("  需要手工或用 tabula-py 提取表格数据")
    print("  安装: pip install tabula-py")
    print("  示例: import tabula; tables = tabula.read_pdf('yearbook.pdf', pages='all')")


def download_phsciencedata():
    """
    公共卫生科学数据中心数据目录和可下载数据
    """
    print("\n" + "="*60)
    print("3. 公共卫生科学数据中心")
    print("="*60)

    print("  网址: https://www.phsciencedata.cn/")
    print("  步骤:")
    print("  1. 注册账号: https://www.phsciencedata.cn/Share/user/login.jsp")
    print("  2. 搜索 '肺结核' 或浏览传染病分类目录")
    print("  3. 下载公开数据集; 超范围数据填写申请表")
    print("  4. 数据申请页: https://www.phsciencedata.cn/Share/edtShare.jsp")
    print("  5. 联系电话: 010-58900412")
    print()
    print("  可获取:")
    print("    - 2004年至今各省年度肺结核发病率(公开)")
    print("    - 月度/县级个案数据(需申请)")


# ============================================================
# 渠道二: 学术论文补充数据
# ============================================================

def download_paper_supplements():
    """
    从已发表论文获取补充数据材料
    """
    print("\n" + "="*60)
    print("4. 学术论文补充数据")
    print("="*60)

    out_dir = CHANNEL2_DIR / "paper_supplements"

    papers = [
        {
            "id": "frontiers_2024_1436515",
            "title": "TB disease burden in China: spatio-temporal clustering (2004-2018, 31省月度)",
            "url": "https://www.frontiersin.org/journals/public-health/articles/10.3389/fpubh.2024.1436515/full",
            "note": "Frontiers OA - 访问论文页面下载 Supplementary Material"
        },
        {
            "id": "frontiers_2026_1866155",
            "title": "Trends and spatial distribution of PTB in China (2004-2025, 31省月度)",
            "url": "https://www.frontiersin.org/journals/public-health/articles/10.3389/fpubh.2026.1866155/full",
            "note": "Frontiers OA - 最新最完整的省级月度数据"
        },
        {
            "id": "sciencedirect_S120197122600322X",
            "title": "Spatiotemporal patterns of PTB in China 2000-2024 (31省月度+预测到2030)",
            "url": "https://www.sciencedirect.com/science/article/pii/S120197122600322X",
            "note": "ScienceDirect - 可能需机构订阅, 覆盖2000-2024"
        },
        {
            "id": "springer_s41182_026_00928",
            "title": "TB incidence and mortality trends in mainland China 2004-2024",
            "url": "https://link.springer.com/article/10.1186/s41182-026-00928-4",
            "note": "Springer OA (Tropical Medicine and Health)"
        },
        {
            "id": "pmc_11149027_jiangsu",
            "title": "PTB spatial-temporal distribution in eastern China (江苏95县 2011-2021)",
            "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11149027/",
            "note": "PMC OA - 县级数据"
        },
        {
            "id": "pmc_9695586_shandong",
            "title": "County-level TB in Shandong 2016-2020 (135,185 cases)",
            "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9695586/",
            "note": "PMC OA - 县级数据"
        },
        {
            "id": "pmc_8297355_hunan",
            "title": "Spatiotemporal patterns of TB in Hunan Province (2013-2018, 县级)",
            "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8297355/",
            "note": "PMC OA - 县级数据"
        },
        {
            "id": "plosone_0290978_mgwr",
            "title": "Influential factors of TB in mainland China (MGWR, 县级面板数据)",
            "url": "https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0290978",
            "note": "PLOS ONE OA - 附件含县级面板数据"
        },
        {
            "id": "pmc_10041483_spatiotemporal",
            "title": "Spatial-temporal dynamic distribution of TB in China 2008-2018 (省级)",
            "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC10041483/",
            "note": "PMC OA - 省级年度面板数据"
        },
        {
            "id": "bmc_s12889_025_24575",
            "title": "Factors affecting TB incidence in mainland Chinese economic regions (省级)",
            "url": "https://link.springer.com/article/10.1186/s12889-025-24575-2",
            "note": "BMC Public Health OA"
        },
    ]

    # Save paper index
    with open(out_dir / "paper_index.json", 'w', encoding='utf-8') as f:
        json.dump(papers, f, ensure_ascii=False, indent=2)

    for p in papers:
        paper_dir = out_dir / p['id']
        paper_dir.mkdir(exist_ok=True)

        print(f"\n  [{p['id']}]")
        print(f"  {p['title']}")
        print(f"  URL: {p['url']}")

        # Try to download the paper page
        try:
            r = requests.get(p['url'], headers=HEADERS, timeout=30)
            with open(paper_dir / "page.html", 'w', encoding='utf-8') as f:
                f.write(r.text)

            # Try to find supplementary material links
            from bs4 import BeautifulSoup
            soup = BeautifulSoup(r.text, 'html.parser')

            supp_links = []
            for a in soup.find_all('a', href=True):
                href = a['href']
                text = a.get_text(strip=True).lower()
                if any(kw in text for kw in ['supplement', 'additional', 'appendix', 'data']):
                    if any(ext in href for ext in ['.xlsx', '.csv', '.xls', '.zip', '.docx', '.pdf']):
                        supp_links.append((a.get_text(strip=True), href))
                elif any(ext in href for ext in ['.xlsx', '.csv', '.xls']):
                    supp_links.append((a.get_text(strip=True), href))

            if supp_links:
                for name, href in supp_links:
                    if not href.startswith('http'):
                        href = p['url'].rsplit('/', 1)[0] + '/' + href
                    ext = href.rsplit('.', 1)[-1][:4]
                    fname = re.sub(r'[^\w]', '_', name)[:50] + '.' + ext
                    download_file(href, paper_dir / fname, name)
                    time.sleep(1)
            else:
                print(f"    未自动发现附件链接, 请手动访问论文页面下载 Supplementary Material")

        except ImportError:
            print("    需要 beautifulsoup4")
        except Exception as e:
            print(f"    {e}")


def download_gbd_data():
    """
    下载GBD全球疾病负担数据(IHME) - 中国省级TB估计
    """
    print("\n" + "="*60)
    print("5. GBD数据 (IHME)")
    print("="*60)

    print("  GBD Results Tool: https://vizhub.healthdata.org/gbd-results/")
    print("  需注册IHME账号后下载")
    print()
    print("  推荐查询参数:")
    print("    - GBD 估计年份: 2021")
    print("    - Location: China (subnational)")
    print("    - Cause: Tuberculosis")
    print("    - Measure: Incidence, Deaths, DALYs")
    print("    - Metric: Rate (per 100,000)")
    print("    - Year: 2000-2021")
    print("    - Age: All ages / age-specific")
    print("    - Sex: Both / Male / Female")
    print()
    print("  下载后保存到: raw/channel2/gbd/")


def download_who_csv():
    """
    下载WHO全球TB数据CSV
    """
    print("\n" + "="*60)
    print("6. WHO全球TB数据库CSV")
    print("="*60)

    out_dir = CHANNEL1_DIR / "who_data"
    out_dir.mkdir(exist_ok=True)

    # WHO TB data download page
    who_files = {
        "TB_burden_countries": "https://extranet.who.int/tme/generateCSV.asp?ds=estimates",
        "TB_notifications": "https://extranet.who.int/tme/generateCSV.asp?ds=notifications",
        "TB_outcomes": "https://extranet.who.int/tme/generateCSV.asp?ds=outcomes",
        "MDR_RR_TB": "https://extranet.who.int/tme/generateCSV.asp?ds=mdr_rr_estimates",
        "TB_budget": "https://extranet.who.int/tme/generateCSV.asp?ds=finance",
    }

    for name, url in who_files.items():
        download_file(url, out_dir / f"{name}.csv", f"WHO {name}")
        time.sleep(2)


# ============================================================
# Main
# ============================================================

def main():
    print("=" * 60)
    print("中国肺结核发病率数据采集脚本")
    print(f"运行时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    os.chdir(BASE_DIR)
    ensure_dirs()

    download_cdc_monthly_reports()
    download_nhc_yearbook()
    download_phsciencedata()
    download_paper_supplements()
    download_gbd_data()
    download_who_csv()

    print("\n" + "=" * 60)
    print("采集完成!")
    print("=" * 60)
    print(f"\n数据保存在: {BASE_DIR}")
    print("\n后续步骤:")
    print("1. 检查 raw/channel2/paper_supplements/ 中的论文附件")
    print("   - 手动从论文页面下载 Supplementary Material (.xlsx/.csv)")
    print("   - 重点: frontiers_2026_1866155 (2004-2025省级月度)")
    print("   - 重点: sciencedirect_S120197122600322X (2000-2024省级月度)")
    print("2. 注册 IHME 账号下载 GBD 数据")
    print("3. 注册公共卫生科学数据中心账号获取详细数据")
    print("4. 用 tabula-py 从年鉴PDF提取表格")
    print("5. 运行 parse_and_merge.py 整合所有数据")


if __name__ == "__main__":
    main()
