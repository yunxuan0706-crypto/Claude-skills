# 中国肺结核发病率数据清单

> 采集日期: 2026-08-30
> 目标: 2000-2025年 省/市/县级 月度及以上 肺结核发病率数据

## 已下��数据 (processed/)

| 文件 | 覆盖范围 | 空间粒度 | 时间粒度 | 来源 |
|------|----------|----------|----------|------|
| `01_china_national_who_annual_estimates_2000_2023.csv` | 2000-2023 | 全国 | 年度 | WHO/TidyTuesday |
| `02_china_national_who_gtb_annual_notifications_2000_2024.csv` | 2000-2024 | 全国 | 年度 | WHO GTB Report 2025 |
| `03_china_national_who_gtb_monthly_ptb_2020_2022.csv` | 2020-2022 | 全国 | 月度 | WHO GTB Report 2025 |
| `04_china_national_who_gtb_full_notifications_1980_2024.csv` | 1980-2024 | 全国 | 年度 | WHO GTB Report 2025 |
| `05_liaoning_monthly_tb_incidence_2005_2017.csv` | 2005-2017 | 辽宁省 | 月度 | EnbinYang/tb_prediction |

## 已下载原始数据 (raw/)

### channel1/ (公开平台)

- `github-infectious-disease/` - 23种乙类传染病月度数据 (2004.04-2020.09)
  - 全国级，含肺结核，Excel格式 (Figure 1-5.xlsx)
  - 数据是分析图表形式，非按病种原始列表
- `tb_prediction_files/` - 辽宁省TB预测研究数据
  - 月度发病率 + 气象/经济协变量 (2005-2017)
- `who-tb-burden/` - WHO TB负担估计值 (TidyTuesday整理)
  - 中国2000-2023年度估计发病率/死亡率
- `China-CDC-raw-data/` - 中国CDC原始数据(登革热/狂犬病等)
  - **不含肺结核** (仅含6种其他疾病)

### channel2/ (论文/国际组织)

- `gtbreport2025/` - WHO全球结核病报告2025完整代码和数据
  - `data/gtb/snapshot_2025-07-30/` 含CSV和RDA格式核心数据
  - annual_notifs.csv: 各国年度通报数 (2000-2024)
  - sub_annual_notifs.csv: 各国月度/季度通报 (中国2020-2022)
  - tb.rda: 完整通报历史 (1980-2024)

## 未能下载 (需本地运行采集脚本)

### 受网络限制无法访问的数据源

| 数据源 | 内容 | 下载方式 |
|--------|------|----------|
| 国家疾控局月报 | 2004-2025全国月度TB发病/死亡 | 运行 `scripts/download_all_sources.py` |
| 中国CDC结核中心月报 | TB专项月报 | 同上 |
| 卫生健康统计年鉴PDF | 31省年度发病率 | 同上 + `scripts/parse_yearbook_pdf.py` |
| Frontiers论文附件 | 31省月度2004-2025 | 同上，手动下载Supplementary |
| ScienceDirect论文附件 | 31省月度2000-2024 | 同上 |
| PMC论文附件 | 江苏/山东/湖南县级数据 | 同上 |
| PLOS ONE附件 | 全国县级面板数据 | 同上 |
| WHO TB数据库CSV | 全球TB通报/估计/治疗结局 | 同上 |
| GBD Results Tool | 中国34省级1990-2023估计 | 注册IHME账号后下载 |
| 公共卫生科学数据中心 | 省/县级详细数据 | 注册后检索下载或申请 |

## 本地运行步骤

```bash
# 1. 安装依赖
pip install requests beautifulsoup4 pandas openpyxl lxml tabula-py

# 2. 运行主采集脚本
cd china-tb-data
python scripts/download_all_sources.py

# 3. 解析CDC月报(采集完成后)
python scripts/parse_cdc_reports.py

# 4. 解析年鉴PDF(需要Java)
python scripts/parse_yearbook_pdf.py

# 5. 手动下载论文附件
# 重点论文(省级月度最完整):
# - https://www.frontiersin.org/articles/10.3389/fpubh.2026.1866155 (2004-2025)
# - https://www.sciencedirect.com/science/article/pii/S120197122600322X (2000-2024)
```

## 数据使用注意

1. **口径差异**: 2004年前为被动报卡体系，2004年起为网络直报，两期数据完整性有系统性差异
2. **概念区分**: 报告发病率(reported incidence) ≠ 调查患病率(prevalence) ≠ WHO估计发病率(estimated incidence)
3. **WHO调整系数**: WHO月度数据乘以了0.7019的调整因子，基于2017-2019年��告数与最终通报数的历史关系
4. **引用**: 使用数据时请引用原始数据来源(WHO/中国CDC/相应论文)
