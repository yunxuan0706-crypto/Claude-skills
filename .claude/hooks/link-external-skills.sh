#!/bin/bash
set -euo pipefail

if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

REPO_ROOT="${CLAUDE_PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
EXTERNAL="$REPO_ROOT/external"
DEST="$HOME/.claude/skills"

cd "$REPO_ROOT"
git submodule update --init --recursive --depth 1 >/dev/null 2>&1 || git submodule update --init --recursive

mkdir -p "$DEST"

link() {
  local dest_name="$1"
  local src="$2"
  if [ -f "$src/SKILL.md" ]; then
    ln -sfn "$src" "$DEST/$dest_name"
  else
    echo "link-external-skills: skipping $dest_name, no SKILL.md at $src" >&2
  fi
}

# 12 academic / writing skills
link "scipilot-figure-skill"     "$EXTERNAL/scipilot-figure-skill"
link "cite-verify"               "$EXTERNAL/claude-skill-citation-checker"
link "stats-sanity"              "$EXTERNAL/claude-statistical-analysis-skill"
link "survey-builder"            "$EXTERNAL/literature-survey-skill"
link "paper-spine"               "$EXTERNAL/PaperSpine/dist/claude/skills/paper-spine"
link "latex-writer"              "$EXTERNAL/latex-precision-skill/.claude/skills/latex-precision-skill"
link "paper-rag"                 "$EXTERNAL/paper-rag/src/paper_rag/data"
link "academic-research"         "$EXTERNAL/academic-research-skills/academic-pipeline"
link "academic-deep-research"    "$EXTERNAL/academic-research-skills/deep-research"
link "academic-paper"            "$EXTERNAL/academic-research-skills/academic-paper"
link "academic-paper-reviewer"   "$EXTERNAL/academic-research-skills/academic-paper-reviewer"
link "ai-research-skill"         "$EXTERNAL/AI-Research-SKILLs/0-autoresearch-skill"

# 42 AI engineering skills from AI-Research-SKILLs, one per category subdirectory
AI_RESEARCH_SKILLS_CATEGORIES=(
  "01-model-architecture"
  "03-fine-tuning"
  "06-post-training"
  "08-distributed-training"
  "15-rag"
  "18-multimodal"
  "12-inference-serving"
)

for category in "${AI_RESEARCH_SKILLS_CATEGORIES[@]}"; do
  category_dir="$EXTERNAL/AI-Research-SKILLs/$category"
  [ -d "$category_dir" ] || continue
  for skill_dir in "$category_dir"/*/; do
    [ -d "$skill_dir" ] || continue
    skill_name="$(basename "$skill_dir")"
    link "$skill_name" "${skill_dir%/}"
  done
done

echo "link-external-skills: linked skills into $DEST"
