"""Generate BUILD_SHEET.md: every shape's exact geometry in centimetres, so the
figure can be rebuilt by hand in PowerPoint's Format Shape pane."""
import json, sys

spec = json.load(open(sys.argv[1] if len(sys.argv) > 1 else 'fig1_ppt_spec.json'))
W, H = spec['w_mm'], spec['h_mm']
cm = lambda v: round(v / 10.0, 2)

BASE = {'#9DC3E6': ('S', 'susceptible'), '#DE6B63': ('I', 'infected'),
        '#43474A': ('R', 'recovered'),   '#FFFFFF': ('U', 'test node u')}
SIZE = {1.55: '', 2.0: '-lg', 1.45: '-m', 1.1: '-sm'}

def node_style(o):
    """Unique style name for a node, from fill / radius / fade."""
    letter, desc = BASE.get(o['fill'], ('?', ''))
    nm = 'N-' + letter + SIZE.get(round(o['r'], 2), '-r%.2f' % o['r'])
    if o['alpha'] < 0.99: nm += '-fade'
    return nm, desc

def he_style(o):
    if o['edge'] == '#B8BCC0': return 'HE-cloud'
    if o['edge_alpha'] < 0.9:  return 'HE-mute'
    return 'HE-b' if o['dashed'] else 'HE-a'

def fade_hex(c, a):
    n = int(c.lstrip('#'), 16)
    mix = lambda k: round(((n >> k) & 255) * a + 255 * (1 - a))
    return '#%02X%02X%02X' % (mix(16), mix(8), mix(0))
PANEL = [('a', 5, 95, 0, 69), ('b', 95, 186, 0, 69),
         ('c', 5, 95, 69, 131), ('d', 95, 186, 69, 131),
         ('legend', 0, 190, 131, 145)]

def panel_of(x, y):
    for name, x0, x1, y0, y1 in PANEL:
        if x0 <= x < x1 and y0 <= y < y1:
            return name
    return 'legend'

def centre(o):
    k = o['kind']
    if k == 'ellipse':  return o['cx'], o['cy']
    if k in ('node', 'cross', 'text'): return o['x'], o['y']
    if k == 'roundbox':  return o['x'] + o['w'] / 2, o['y'] + o['h'] / 2
    return (o['p0'][0] + o['p1'][0]) / 2, (o['p0'][1] + o['p1'][1]) / 2

groups = {n: {'ellipse': [], 'node': [], 'line': [], 'text': [], 'box': []}
          for n, *_ in PANEL}
for o in spec['shapes']:
    g = groups[panel_of(*centre(o))]
    if o['kind'] == 'ellipse':                       g['ellipse'].append(o)
    elif o['kind'] == 'node':                        g['node'].append(o)
    elif o['kind'] in ('arrow', 'rule', 'cross'):    g['line'].append(o)
    elif o['kind'] == 'roundbox':                    g['box'].append(o)
    else:                                            g['text'].append(o)

def runs_txt(rs):
    """Group consecutive runs by baseline so (a) reads as one superscript."""
    out, buf, state = [], '', None
    def flush():
        nonlocal buf
        if buf:
            out.append({'sup': '^{%s}', 'sub': '_{%s}'}.get(state, '%s') % buf)
        buf = ''
    for r in rs:
        st = 'sup' if r['sup'] else 'sub' if r['sub'] else None
        if st != state:
            flush(); state = st
        buf += r['t']
    flush()
    return ''.join(out)

L = []
w = L.append
w('# Fig. 1 — PowerPoint 施工表\n')
w('本表对应 `fig1_ppt_reproducible.pptx` / `fig1_ppt_190mm.png`（照着画的目标图）。\n')
w('**最快的做法**：直接打开 `fig1_ppt_reproducible.pptx`，把里面的形状复制到你自己的'
  '幻灯片里再改。本表是给需要从零手画、或要核对精确数值时用的。\n')
w('## 0. 画布与通用设置\n')
w('| 项目 | 值 |\n|---|---|')
w('| 设计 → 幻灯片大小 → 自定义 | **%.2f cm × %.2f cm** |' % (cm(W), cm(H)))
w('| 字体（全图统一） | **Arial** |')
w('| 单位 | 厘米。表中「X / Y」= 形状外框左上角，填在「设置形状格式 → 大小与属性 → 位置」 |')
w('| 超边、节点 | 插入 → 形状 → **椭圆**（按住 Shift 得正圆） |')
w('| 箭头、分隔线 | 插入 → 形状 → **直线**；箭头在「线条 → 箭头末端类型 → 箭头」 |')
w('| 公式框 | 插入 → 形状 → **圆角矩形** |')
w('| 颜色 | 填充/线条 → 其他颜色 → 自定义 → 十六进制 |\n')

w('## 1. 样式表（先定义，后面的表只引用名字）\n')
w('### 节点（正圆）\n')
w('| 样式 | 含义 | 直径 cm | 填充 | 线条色 | 线宽 pt |\n|---|---|---|---|---|---|')
seen = {}
for o in spec['shapes']:
    if o['kind'] == 'node':
        seen.setdefault(node_style(o)[0], o)
for nm in sorted(seen):
    o = seen[nm]; desc = node_style(o)[1]
    a = o['alpha']
    w('| `%s` | %s%s | %.2f | `%s` | `%s` | %.2f |'
      % (nm, desc, '' if a > 0.99 else '（淡化：填充 `%s`，线条 `%s`）'
         % (fade_hex(o['fill'], a), fade_hex(o['edge'], a)),
         cm(2 * o['r']), o['fill'], o['edge'], o['lw']))
w('')
w('### 超边（椭圆，只描边 + 极淡填充）\n')
w('| 样式 | 用于 | 填充 | 填充透明度 | 线条色 | 线宽 pt | 虚线 |\n'
  '|---|---|---|---|---|---|---|')
DESC = {'HE-a': 'layer $a$ 超边', 'HE-b': 'layer $b$ 超边',
        'HE-mute': '已到达超边（排除，整体淡化）', 'HE-cloud': 'rest of network'}
hs = {}
for o in spec['shapes']:
    if o['kind'] == 'ellipse':
        hs.setdefault(he_style(o), o)
for nm in ('HE-a', 'HE-b', 'HE-mute', 'HE-cloud'):
    if nm not in hs: continue
    o = hs[nm]
    w('| `%s` | %s | `%s` | %d%% | `%s` | %.2f | %s |'
      % (nm, DESC[nm], o['fill'], round((1 - o['alpha']) * 100),
         fade_hex(o['edge'], o['edge_alpha']), o['lw'],
         '**划线**' if o['dashed'] else '实线'))
w('')
w('### 线条与文字\n')
w('| 样式 | 用途 | 颜色 | 线宽 pt | 箭头 |\n|---|---|---|---|---|')
w('| `AR-T` | 传播箭头 | `#C0392B` | 0.80–0.85 | 末端箭头 |')
w('| `AR-N` | 闭合环 / 级联箭头 | `#4A4A4A` | 0.80 | 末端箭头 |')
w('| `AR-C` | $h_a$ 耦合箭头 | `#3D3D3D` | 0.85 | 末端箭头 |')
w('| `RULE` | 面板分隔线 | `#DCDCDC` | 0.60 | 无 |')
w('| `X` | 排除标记（两条交叉线） | `#606060` | 1.40 | 无 |\n')
w('| 文字样式 | 字号 pt | 粗体 |\n|---|---|---|')
w('| 面板标号 `(a)` | 9.0 | 是 |')
w('| 面板小标题 | 7.5 | 否 |')
w('| 符号 / 公式 | 7.8 | 否（$K_{ab}$ 那格 8.6 粗体） |')
w('| 行标签、图例 | 7.0 | 否 |\n')

TITLE = {'a': '面板 (a) Multiplex hypergraph', 'b': '面板 (b) Edge-based closure',
         'c': '面板 (c) Within-hyperedge factor C',
         'd': '面板 (d) Excess hyperdegree B', 'legend': '图例与分隔线'}
for n, *_ in PANEL:
    g = groups[n]
    w('## %s\n' % TITLE[n])
    if g['ellipse']:
        w('**超边（椭圆）** — 宽/高是外框尺寸；旋转填在「大小与属性 → 旋转」\n')
        w('| # | X cm | Y cm | 宽 cm | 高 cm | 旋转° | 样式 |\n|---|---|---|---|---|---|---|')
        for i, o in enumerate(g['ellipse'], 1):
            st = he_style(o)
            w('| %d | %.2f | %.2f | %.2f | %.2f | %.1f | `%s` |'
              % (i, cm(o['cx'] - o['a']), cm(o['cy'] - o['b']),
                 cm(2 * o['a']), cm(2 * o['b']), o['ang'], st))
        w('')
    if g['node']:
        w('**节点** — X/Y 是外框左上角\n')
        w('| # | X cm | Y cm | 样式 |\n|---|---|---|---|')
        for i, o in enumerate(g['node'], 1):
            w('| %d | %.2f | %.2f | `%s` |'
              % (i, cm(o['x'] - o['r']), cm(o['y'] - o['r']), node_style(o)[0]))
        w('')
    if g['line']:
        w('**直线 / 箭头** — 起点 → 终点；PowerPoint 里画好后按 X/Y/宽/高 对位\n')
        w('| # | 起点 X,Y cm | 终点 X,Y cm | X cm | Y cm | 宽 cm | 高 cm | 样式 |\n'
          '|---|---|---|---|---|---|---|---|')
        i = 0
        for o in g['line']:
            segs = ([((o['x'] - o['s'], o['y'] - o['s']), (o['x'] + o['s'], o['y'] + o['s'])),
                     ((o['x'] - o['s'], o['y'] + o['s']), (o['x'] + o['s'], o['y'] - o['s']))]
                    if o['kind'] == 'cross' else [(o['p0'], o['p1'])])
            st = ('X' if o['kind'] == 'cross' else 'RULE' if o['kind'] == 'rule'
                  else 'AR-T' if o['color'] == '#C0392B'
                  else 'AR-C' if o['color'] == '#3D3D3D' else 'AR-N')
            for p0, p1 in segs:
                i += 1
                w('| %d | %.2f, %.2f | %.2f, %.2f | %.2f | %.2f | %.2f | %.2f | `%s` |'
                  % (i, cm(p0[0]), cm(p0[1]), cm(p1[0]), cm(p1[1]),
                     cm(min(p0[0], p1[0])), cm(min(p0[1], p1[1])),
                     cm(abs(p1[0] - p0[0])), cm(abs(p1[1] - p0[1])), st))
        w('')
    if g['box']:
        w('**圆角矩形**\n')
        w('| X cm | Y cm | 宽 cm | 高 cm | 填充 | 线条 |\n|---|---|---|---|---|---|')
        for o in g['box']:
            w('| %.2f | %.2f | %.2f | %.2f | `%s` | `%s` 0.70 pt |'
              % (cm(o['x']), cm(o['y']), cm(o['w']), cm(o['h']), o['fill'], o['edge']))
        w('')
    if g['text']:
        w('**文本框** — 内边距设 0；`^{}` = 上标，`_{}` = 下标，斜体按数学惯例\n')
        w('| # | X cm | Y cm | 对齐 | 字号 pt | 内容 |\n|---|---|---|---|---|---|')
        for i, o in enumerate(g['text'], 1):
            w('| %d | %.2f | %.2f | %s | %.1f | `%s` |'
              % (i, cm(o['x']), cm(o['y']),
                 {'left': '左', 'center': '中', 'right': '右'}[o['ha']],
                 o['size'], runs_txt(o['runs'])))
        w('')

open('BUILD_SHEET.md', 'w').write('\n'.join(L))
print('wrote BUILD_SHEET.md  (%d lines)' % len(L))
