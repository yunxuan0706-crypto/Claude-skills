# Claude Skills

个人 Claude Code 技能集合。

## 本仓库内置技能

以下技能以 `SKILL.md` 形式收录在本仓库中,可直接在 Claude Code 中使用:

| 技能 | 说明 |
| --- | --- |
| [`coding`](coding/) | 编码任务通用工作流 |
| [`data-analysis`](data-analysis/) | 数据分析工作流 |
| [`latex-manuscript`](latex-manuscript/) | LaTeX 学术论文撰写、排版与排错 |
| [`literature-review`](literature-review/) | 文献综述与相关工作梳理 |
| [`peer-review-response`](peer-review-response/) | 同行评审意见回复 |
| [`research`](research/) | 多源信息研究与综合 |
| [`scientific-writing`](scientific-writing/) | 科研论文写作 |

## 外部技能(External Skills,git submodule)

以下是 GitHub 上与学术科研工作流相关的第三方 Claude Skill 项目,以 **git submodule** 的形式接入 `external/` 目录 —— 不复制/vendor 其源码,只保留指向对应上游仓库指定 commit 的指针,更新时可随时对齐上游最新版本。

| 类别 | 子模块路径 | 上游仓库 | 说明 |
| --- | --- | --- | --- |
| 科研图表 | [`external/scipilot-figure-skill`](external/scipilot-figure-skill) | [Haojae/scipilot-figure-skill](https://github.com/Haojae/scipilot-figure-skill) | 发表级科研图表 copilot,先诊断数据、推荐图表类型,再产出符合 Nature/Science/IEEE 等期刊规范的图 |
| 自动化科研工具集 | [`external/awesome-auto-research-tools`](external/awesome-auto-research-tools) | [handsome-rich/Awesome-Auto-Research-Tools](https://github.com/handsome-rich/Awesome-Auto-Research-Tools) | 自动化科研工具精选列表,涵盖文献检索、读文献、实验管理、代码生成 |
| 学术研究全流程 | [`external/academic-research-skills`](external/academic-research-skills) | [Imbad0202/academic-research-skills](https://github.com/imbad0202/academic-research-skills) | 覆盖 研究 → 写作 → 审阅 → 修订 → 定稿 全流程的学术研究技能套件 |
| 引用核验 | [`external/claude-skill-citation-checker`](external/claude-skill-citation-checker) | [PHY041/claude-skill-citation-checker](https://github.com/PHY041/claude-skill-citation-checker) | 将 `.bib` 引用与 CrossRef / Semantic Scholar / OpenAlex 交叉核验,检测幻觉引用 |
| 论文 RAG 检索 | [`external/paper-rag`](external/paper-rag) | [LucasJLBraz/paper-rag](https://github.com/LucasJLBraz/paper-rag) | 本地 RAG 系统,对一整个文件夹的 PDF 论文做混合检索(向量 + BM25),无需每次重读全文 |
| 论文主线梳理 | [`external/PaperSpine`](external/PaperSpine) | [WUBING2023/PaperSpine](https://github.com/WUBING2023/PaperSpine) | 从强论文中提炼核心论证主线,辅助论文改写与证据驱动的修订 |
| LaTeX 写作 | [`external/latex-precision-skill`](external/latex-precision-skill) | [MagicMonkey-XK/latex-precision-skill](https://github.com/MagicMonkey-XK/latex-precision-skill) | 逐行精修的 LaTeX 技能,支持公式、表格、交叉引用、参考文献与中文排版 |
| 统计结果自检 | [`external/claude-statistical-analysis-skill`](external/claude-statistical-analysis-skill) | [TerryFYL/claude-statistical-analysis-skill](https://github.com/TerryFYL/claude-statistical-analysis-skill) | 分析前先诊断数据质量与统计假设是否成立,自动选择合适方法,输出 APA 格式结果 |
| 文献综述构建 | [`external/literature-survey-skill`](external/literature-survey-skill) | [SNL-UCSB/literature-survey-skill](https://github.com/SNL-UCSB/literature-survey-skill) | 把一堆论文整理为对该领域"已知/未知/可推进方向"的系统性综述 |
| AI 科研技能库 | [`external/AI-Research-SKILLs`](external/AI-Research-SKILLs) | [Orchestra-Research/AI-Research-SKILLs](https://github.com/Orchestra-Research/AI-Research-SKILLs) | 覆盖模型架构、训练、评测、RAG、Agent 等方向的大型开源 AI 科研技能库 |

> 以上外部项目均为第三方维护,子模块只固定了当前接入时的 commit,使用前请查看各自仓库的 License 与更新状态。

### 如何使用这些外部技能

**1. clone 本仓库时带上子模块内容**

```bash
git clone --recurse-submodules <本仓库地址>
# 如果已经 clone 过、忘了带 --recurse-submodules:
git submodule update --init --recursive
```

**2. 让 Claude Code 能发现这些技能**

Claude Code 会扫描 `~/.claude/skills/`(全局)或项目内 `.claude/skills/`(项目级)下含有 `SKILL.md` 的目录。把需要的子模块目录软链接或复制过去即可,例如:

```bash
# 全局可用(推荐,任何项目里都能触发)
ln -s "$(pwd)/external/paper-rag" ~/.claude/skills/paper-rag
ln -s "$(pwd)/external/claude-skill-citation-checker" ~/.claude/skills/cite-verify
# 其余 8 个同理,把 external/<name> 软链到 ~/.claude/skills/<你喜欢的名字>
```

**3. 调用方式**

- 多数技能是**语义触发**:直接在 Claude Code 对话里描述意图(比如"帮我核对这篇论文的引用是否准确"、"根据这个 CSV 画一张适合投 IEEE 的图"),Claude 会自动匹配到对应技能,不需要记命令。
- 部分技能提供**显式斜杠命令**(如 `/verify-citations`、`/survey` 等),具体命令名以各自仓库根目录的 `SKILL.md`(或 `.claude/commands/`)为准,进入对应 `external/<name>` 目录查看即可。
- `Awesome-Auto-Research-Tools` 和 `AI-Research-SKILLs` 体量较大、覆盖面广,建议先打开其 `README.md` 挑选实际需要的子技能再软链,而不是整体接入,避免技能过多导致触发冲突。

**4. 更新到上游最新版本**

```bash
git submodule update --remote external/<name>
git add external/<name>
git commit -m "Bump <name> submodule"
```
