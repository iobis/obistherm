<script lang="ts">
  import { db, runQuery } from '$lib/db.svelte';

  let {
    studyId = null as string | null,
    onStudyChange = (_id: string) => {},
  }: {
    studyId: string | null;
    onStudyChange: (id: string) => void;
  } = $props();

  // ── Study definitions ─────────────────────────────────────────────────────
  const STUDIES = [
    {
      id: 'fiddler-warming',
      title: 'Fiddler crabs & warming',
      subtitle: 'Temperature experienced by Leptuca thayeri and Minuca rapax (2000–2024)',
      description: 'Two intertidal fiddler crab species from the western Atlantic illustrate how even a modest warming signal (~0.5°C over 25 years) becomes visible in occurrence-matched temperature records.',
      sql: `SELECT year, species, ROUND(AVG(coraltempSST), 3) AS mean_sst, COUNT(*) AS records
FROM obistherm
WHERE species IN ('Leptuca thayeri', 'Minuca rapax')
GROUP BY year, species
ORDER BY year, species`,
    },
    {
      id: 'sst-agreement',
      title: 'SST product agreement',
      subtitle: 'CoralTemp vs MUR for all observations, grouped by family',
      description: 'Four independent SST products — GLORYS, CoralTemp, MUR, and OSTIA — track the same physical signal but differ in resolution and methodology. The synthetic data encodes small random offsets between products to illustrate real-world biases.',
      sql: `SELECT ROUND(coraltempSST, 1) AS coraltemp, ROUND(murSST, 1) AS mur, family
FROM obistherm
WHERE coraltempSST IS NOT NULL AND murSST IS NOT NULL
ORDER BY coraltemp`,
    },
    {
      id: 'seasonality',
      title: 'Seasonal temperature patterns',
      subtitle: 'Monthly mean SST by habitat group',
      description: 'Tropical reef species experience warm, stable temperatures year-round. Cold-water species show pronounced summer peaks. Subtropical coastal species fall between these extremes — a classic biodiversity-temperature interaction.',
      sql: `SELECT month,
  CASE
    WHEN family IN ('Acanthuridae', 'Acroporidae') THEN 'Tropical reef'
    WHEN family IN ('Gadidae', 'Scombridae') THEN 'Cold-temperate ocean'
    ELSE 'Subtropical coast'
  END AS habitat,
  ROUND(AVG(coraltempSST), 2) AS mean_sst
FROM obistherm
WHERE coraltempSST IS NOT NULL
GROUP BY month, habitat
ORDER BY month, habitat`,
    },
  ];

  // ── Query state ───────────────────────────────────────────────────────────
  type QueryResult = Record<string, unknown>[];
  let results = $state<Record<string, QueryResult>>({});
  let loading = $state<Record<string, boolean>>({});
  let errors  = $state<Record<string, string>>({});

  async function runStudy(id: string) {
    const study = STUDIES.find(s => s.id === id);
    if (!study) return;
    loading[id] = true;
    errors[id]  = '';
    try {
      results[id] = await runQuery(study.sql);
    } catch (e) {
      errors[id] = String(e);
    } finally {
      loading[id] = false;
    }
  }

  // Auto-run when studyId changes and db is ready
  $effect(() => {
    if (studyId && db.status === 'ready' && !results[studyId] && !loading[studyId]) {
      runStudy(studyId);
    }
  });

  // ── Chart dimensions ──────────────────────────────────────────────────────
  const LINE_PAD  = { top: 24, right: 130, bottom: 48, left: 58 };
  const LINE_W    = 600;
  const LINE_H    = 300;
  const LINE_PW   = LINE_W - LINE_PAD.left - LINE_PAD.right; // 412
  const LINE_PH   = LINE_H - LINE_PAD.top  - LINE_PAD.bottom; // 228

  const SCAT_PAD  = { top: 24, right: 24, bottom: 80, left: 58 };
  const SCAT_W    = 600;
  const SCAT_H    = 380;
  const SCAT_PW   = SCAT_W - SCAT_PAD.left - SCAT_PAD.right;
  const SCAT_PH   = SCAT_H - SCAT_PAD.top  - SCAT_PAD.bottom;

  // ── Fiddler crabs chart ───────────────────────────────────────────────────
  const FIDDLER_COLORS: Record<string, string> = {
    'Leptuca thayeri': '#f68512',
    'Minuca rapax':    '#11b5ae',
  };

  const fiddlerData = $derived(() => {
    const rows = results['fiddler-warming'] ?? [];
    if (!rows.length) return null;

    const bySpecies: Record<string, {year: number; mean_sst: number}[]> = {};
    for (const row of rows) {
      const sp  = String(row.species);
      const yr  = Number(row.year);
      const sst = Number(row.mean_sst);
      if (!bySpecies[sp]) bySpecies[sp] = [];
      bySpecies[sp].push({ year: yr, mean_sst: sst });
    }

    const allYears = rows.map(r => Number(r.year));
    const allSSTs  = rows.map(r => Number(r.mean_sst));
    const minY = Math.min(...allYears), maxY = Math.max(...allYears);
    const minSST = Math.floor(Math.min(...allSSTs) - 0.5);
    const maxSST = Math.ceil(Math.max(...allSSTs)  + 0.5);

    function xScale(y: number) { return LINE_PAD.left + (y - minY) / (maxY - minY || 1) * LINE_PW; }
    function yScale(s: number) { return LINE_PAD.top  + (1 - (s - minSST) / (maxSST - minSST || 1)) * LINE_PH; }

    const species = Object.keys(bySpecies);
    const lines = species.map(sp => {
      const pts = bySpecies[sp].sort((a, b) => a.year - b.year);
      const points = pts.map(p => `${xScale(p.year)},${yScale(p.mean_sst)}`).join(' ');
      const dots   = pts.map(p => ({ cx: xScale(p.year), cy: yScale(p.mean_sst) }));
      return { sp, color: FIDDLER_COLORS[sp] ?? '#94a3b8', points, dots };
    });

    // Axes
    const xTicks = Array.from({ length: maxY - minY + 1 }, (_, i) => minY + i)
      .filter(y => y % 5 === 0);
    const yTicks = Array.from({ length: Math.round((maxSST - minSST) / 0.5) + 1 }, (_, i) =>
      Math.round((minSST + i * 0.5) * 10) / 10
    );

    return { lines, xTicks, yTicks, xScale, yScale, minSST, maxSST, minY, maxY, species };
  });

  // ── SST agreement chart ───────────────────────────────────────────────────
  const FAMILY_COLORS: Record<string, string> = {
    'Ocypodidae':   '#f68512',
    'Acanthuridae': '#11b5ae',
    'Gadidae':      '#1473e6',
    'Delphinidae':  '#9b59b6',
    'Acroporidae':  '#e84393',
    'Scombridae':   '#16a34a',
    'Gadiformes':   '#dc2626',
  };
  const FAMILY_PALETTE = ['#f68512','#11b5ae','#1473e6','#9b59b6','#e84393','#16a34a','#dc2626','#64748b'];

  const sstAgreementData = $derived(() => {
    const rows = results['sst-agreement'] ?? [];
    if (!rows.length) return null;

    const allCT  = rows.map(r => Number(r.coraltemp));
    const allMUR = rows.map(r => Number(r.mur));
    const allVals = [...allCT, ...allMUR];
    const minV = Math.floor(Math.min(...allVals) - 1);
    const maxV = Math.ceil(Math.max(...allVals)  + 1);

    function xScale(v: number) { return SCAT_PAD.left + (v - minV) / (maxV - minV) * SCAT_PW; }
    function yScale(v: number) { return SCAT_PAD.top  + (1 - (v - minV) / (maxV - minV)) * SCAT_PH; }

    // Assign colors by family
    const families = [...new Set(rows.map(r => String(r.family)))];
    const familyColor: Record<string, string> = {};
    families.forEach((f, i) => { familyColor[f] = FAMILY_PALETTE[i % FAMILY_PALETTE.length]; });

    const points = rows.map(r => ({
      cx: xScale(Number(r.coraltemp)),
      cy: yScale(Number(r.mur)),
      color: familyColor[String(r.family)] ?? '#94a3b8',
      family: String(r.family),
    }));

    // 1:1 reference line
    const refX1 = xScale(minV), refY1 = yScale(minV);
    const refX2 = xScale(maxV), refY2 = yScale(maxV);

    const ticks = Array.from({ length: Math.round((maxV - minV) / 2) + 1 }, (_, i) => minV + i * 2);

    return { points, families, familyColor, refX1, refY1, refX2, refY2, xScale, yScale, ticks, minV, maxV };
  });

  // ── Seasonality chart ─────────────────────────────────────────────────────
  const HABITAT_COLORS: Record<string, string> = {
    'Tropical reef':        '#e84393',
    'Cold-temperate ocean': '#1473e6',
    'Subtropical coast':    '#f68512',
  };
  const MONTH_NAMES = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'];

  const seasonalityData = $derived(() => {
    const rows = results['seasonality'] ?? [];
    if (!rows.length) return null;

    const byHabitat: Record<string, {month: number; mean_sst: number}[]> = {};
    for (const row of rows) {
      const h   = String(row.habitat);
      const mo  = Number(row.month);
      const sst = Number(row.mean_sst);
      if (!byHabitat[h]) byHabitat[h] = [];
      byHabitat[h].push({ month: mo, mean_sst: sst });
    }

    const allSSTs = rows.map(r => Number(r.mean_sst));
    const minSST  = Math.floor(Math.min(...allSSTs) - 1);
    const maxSST  = Math.ceil(Math.max(...allSSTs)  + 1);

    function xScale(m: number) { return LINE_PAD.left + (m - 1) / 11 * LINE_PW; }
    function yScale(s: number) { return LINE_PAD.top  + (1 - (s - minSST) / (maxSST - minSST || 1)) * LINE_PH; }

    const habitats = Object.keys(byHabitat);
    const lines = habitats.map(h => {
      const pts    = byHabitat[h].sort((a, b) => a.month - b.month);
      const points = pts.map(p => `${xScale(p.month)},${yScale(p.mean_sst)}`).join(' ');
      const dots   = pts.map(p => ({ cx: xScale(p.month), cy: yScale(p.mean_sst) }));
      return { habitat: h, color: HABITAT_COLORS[h] ?? '#94a3b8', points, dots };
    });

    const yTicks = Array.from({ length: Math.round((maxSST - minSST) / 2) + 1 }, (_, i) => minSST + i * 2);

    return { lines, xScale, yScale, yTicks, minSST, maxSST };
  });

  // ── Grid lines ────────────────────────────────────────────────────────────
  function hGridLines(yTicks: number[], yScale: (v: number) => number, width: number, padLeft: number) {
    return yTicks.map(t => {
      const y = yScale(t);
      return { y, x1: padLeft, x2: padLeft + width };
    });
  }

  const activeStudy = $derived(() => STUDIES.find(s => s.id === studyId) ?? null);
</script>

<div class="panel">

  <!-- Study selector buttons -->
  <div class="study-selector">
    {#each STUDIES as s}
      <button
        class="study-btn"
        class:active={studyId === s.id}
        onclick={() => onStudyChange(s.id)}
      >
        {s.title}
      </button>
    {/each}
  </div>

  {#if !studyId}
    <div class="empty">
      <div class="empty-icon">◈</div>
      <p>Select a case study above to explore the data.</p>
    </div>

  {:else if loading[studyId]}
    <div class="empty">
      <div class="spinner"></div>
      <p>Running query…</p>
    </div>

  {:else if errors[studyId]}
    <div class="empty">
      <p class="err">{errors[studyId]}</p>
    </div>

  {:else if activeStudy()}
    {@const study = activeStudy()!}
    <div class="study-body">
      <div class="study-header">
        <div class="study-meta">
          <h2 class="study-title">{study.title}</h2>
          <p class="study-subtitle">{study.subtitle}</p>
          <p class="study-desc">{study.description}</p>
        </div>

        <!-- ── Chart area ── -->
        <div class="chart-area">

          <!-- FIDDLER CRABS LINE CHART -->
          {#if studyId === 'fiddler-warming' && fiddlerData()}
            {@const d = fiddlerData()!}
            <svg viewBox="0 0 {LINE_W} {LINE_H}" class="chart-svg">
              <!-- Grid lines -->
              {#each hGridLines(d.yTicks, d.yScale, LINE_PW, LINE_PAD.left) as gl}
                <line x1={gl.x1} y1={gl.y} x2={gl.x2} y2={gl.y}
                  stroke="#f1f5f9" stroke-width="1" stroke-dasharray="3,3" />
              {/each}

              <!-- Axes -->
              <line x1={LINE_PAD.left} y1={LINE_PAD.top} x2={LINE_PAD.left} y2={LINE_PAD.top + LINE_PH}
                stroke="#e2e8f0" stroke-width="1" />
              <line x1={LINE_PAD.left} y1={LINE_PAD.top + LINE_PH} x2={LINE_PAD.left + LINE_PW} y2={LINE_PAD.top + LINE_PH}
                stroke="#e2e8f0" stroke-width="1" />

              <!-- Y axis ticks + labels -->
              {#each d.yTicks as t}
                {@const y = d.yScale(t)}
                <line x1={LINE_PAD.left - 4} y1={y} x2={LINE_PAD.left} y2={y} stroke="#94a3b8" stroke-width="1" />
                <text x={LINE_PAD.left - 7} y={y + 4} text-anchor="end" font-size="10" fill="#94a3b8"
                  font-family="IBM Plex Mono, monospace">{t.toFixed(1)}</text>
              {/each}

              <!-- X axis ticks + labels -->
              {#each d.xTicks as yr}
                {@const x = d.xScale(yr)}
                <line x1={x} y1={LINE_PAD.top + LINE_PH} x2={x} y2={LINE_PAD.top + LINE_PH + 4}
                  stroke="#94a3b8" stroke-width="1" />
                <text x={x} y={LINE_PAD.top + LINE_PH + 16} text-anchor="middle" font-size="10"
                  fill="#94a3b8" font-family="IBM Plex Mono, monospace">{yr}</text>
              {/each}

              <!-- Y axis label -->
              <text x={14} y={LINE_PAD.top + LINE_PH / 2} text-anchor="middle"
                font-size="10" fill="#64748b" font-family="IBM Plex Mono, monospace"
                transform="rotate(-90, 14, {LINE_PAD.top + LINE_PH / 2})">SST (°C)</text>

              <!-- Lines -->
              {#each d.lines as line}
                <polyline points={line.points} fill="none" stroke={line.color} stroke-width="2" />
                {#each line.dots as dot}
                  <circle cx={dot.cx} cy={dot.cy} r="4" fill={line.color} />
                {/each}
              {/each}

              <!-- Legend -->
              {#each d.species as sp, i}
                <circle cx={LINE_PAD.left + LINE_PW + 14} cy={LINE_PAD.top + 20 + i * 22} r="5"
                  fill={FIDDLER_COLORS[sp] ?? '#94a3b8'} />
                <text x={LINE_PAD.left + LINE_PW + 22} y={LINE_PAD.top + 24 + i * 22}
                  font-size="10" fill="#475569" font-family="IBM Plex Mono, monospace"
                  font-style="italic">{sp}</text>
              {/each}
            </svg>

          <!-- SST AGREEMENT SCATTER PLOT -->
          {:else if studyId === 'sst-agreement' && sstAgreementData()}
            {@const d = sstAgreementData()!}
            <svg viewBox="0 0 {SCAT_W} {SCAT_H}" class="chart-svg">
              <!-- Grid -->
              {#each d.ticks as t}
                {@const y = d.yScale(t)}
                <line x1={SCAT_PAD.left} y1={y} x2={SCAT_PAD.left + SCAT_PW} y2={y}
                  stroke="#f1f5f9" stroke-width="1" stroke-dasharray="3,3" />
              {/each}

              <!-- 1:1 reference line -->
              <line x1={d.refX1} y1={d.refY1} x2={d.refX2} y2={d.refY2}
                stroke="#94a3b8" stroke-width="1.5" stroke-dasharray="6,4" />

              <!-- Axes -->
              <line x1={SCAT_PAD.left} y1={SCAT_PAD.top} x2={SCAT_PAD.left} y2={SCAT_PAD.top + SCAT_PH}
                stroke="#e2e8f0" stroke-width="1" />
              <line x1={SCAT_PAD.left} y1={SCAT_PAD.top + SCAT_PH} x2={SCAT_PAD.left + SCAT_PW} y2={SCAT_PAD.top + SCAT_PH}
                stroke="#e2e8f0" stroke-width="1" />

              <!-- X ticks + labels -->
              {#each d.ticks as t}
                {@const x = d.xScale(t)}
                <line x1={x} y1={SCAT_PAD.top + SCAT_PH} x2={x} y2={SCAT_PAD.top + SCAT_PH + 4}
                  stroke="#94a3b8" stroke-width="1" />
                <text x={x} y={SCAT_PAD.top + SCAT_PH + 16} text-anchor="middle" font-size="10"
                  fill="#94a3b8" font-family="IBM Plex Mono, monospace">{t}</text>
              {/each}

              <!-- Y ticks + labels -->
              {#each d.ticks as t}
                {@const y = d.yScale(t)}
                <line x1={SCAT_PAD.left - 4} y1={y} x2={SCAT_PAD.left} y2={y} stroke="#94a3b8" stroke-width="1" />
                <text x={SCAT_PAD.left - 7} y={y + 4} text-anchor="end" font-size="10"
                  fill="#94a3b8" font-family="IBM Plex Mono, monospace">{t}</text>
              {/each}

              <!-- Axis labels -->
              <text x={SCAT_PAD.left + SCAT_PW / 2} y={SCAT_PAD.top + SCAT_PH + 36}
                text-anchor="middle" font-size="11" fill="#64748b" font-family="IBM Plex Mono, monospace">
                CoralTemp SST (°C)
              </text>
              <text x={14} y={SCAT_PAD.top + SCAT_PH / 2} text-anchor="middle"
                font-size="11" fill="#64748b" font-family="IBM Plex Mono, monospace"
                transform="rotate(-90, 14, {SCAT_PAD.top + SCAT_PH / 2})">MUR SST (°C)</text>

              <!-- Points -->
              {#each d.points as pt}
                <circle cx={pt.cx} cy={pt.cy} r="5" fill={pt.color} opacity="0.65" />
              {/each}

              <!-- Legend (below chart, horizontal) -->
              {#each d.families as fam, i}
                {@const lx = SCAT_PAD.left + (i % 4) * 145}
                {@const ly = SCAT_PAD.top + SCAT_PH + 52 + Math.floor(i / 4) * 18}
                <circle cx={lx + 5} cy={ly} r="5" fill={d.familyColor[fam]} opacity="0.8" />
                <text x={lx + 14} y={ly + 4} font-size="10" fill="#475569"
                  font-family="IBM Plex Mono, monospace">{fam}</text>
              {/each}
            </svg>

          <!-- SEASONALITY LINE CHART -->
          {:else if studyId === 'seasonality' && seasonalityData()}
            {@const d = seasonalityData()!}
            <svg viewBox="0 0 {LINE_W} {LINE_H}" class="chart-svg">
              <!-- Grid -->
              {#each hGridLines(d.yTicks, d.yScale, LINE_PW, LINE_PAD.left) as gl}
                <line x1={gl.x1} y1={gl.y} x2={gl.x2} y2={gl.y}
                  stroke="#f1f5f9" stroke-width="1" stroke-dasharray="3,3" />
              {/each}

              <!-- Axes -->
              <line x1={LINE_PAD.left} y1={LINE_PAD.top} x2={LINE_PAD.left} y2={LINE_PAD.top + LINE_PH}
                stroke="#e2e8f0" stroke-width="1" />
              <line x1={LINE_PAD.left} y1={LINE_PAD.top + LINE_PH} x2={LINE_PAD.left + LINE_PW} y2={LINE_PAD.top + LINE_PH}
                stroke="#e2e8f0" stroke-width="1" />

              <!-- Y ticks -->
              {#each d.yTicks as t}
                {@const y = d.yScale(t)}
                <line x1={LINE_PAD.left - 4} y1={y} x2={LINE_PAD.left} y2={y} stroke="#94a3b8" stroke-width="1" />
                <text x={LINE_PAD.left - 7} y={y + 4} text-anchor="end" font-size="10"
                  fill="#94a3b8" font-family="IBM Plex Mono, monospace">{t}</text>
              {/each}

              <!-- X month labels -->
              {#each MONTH_NAMES as mo, idx}
                {@const x = d.xScale(idx + 1)}
                <text x={x} y={LINE_PAD.top + LINE_PH + 16} text-anchor="middle" font-size="10"
                  fill="#94a3b8" font-family="IBM Plex Mono, monospace">{mo}</text>
              {/each}

              <!-- Y label -->
              <text x={14} y={LINE_PAD.top + LINE_PH / 2} text-anchor="middle"
                font-size="10" fill="#64748b" font-family="IBM Plex Mono, monospace"
                transform="rotate(-90, 14, {LINE_PAD.top + LINE_PH / 2})">SST (°C)</text>

              <!-- Lines -->
              {#each d.lines as line}
                <polyline points={line.points} fill="none" stroke={line.color} stroke-width="2" />
                {#each line.dots as dot}
                  <circle cx={dot.cx} cy={dot.cy} r="4" fill={line.color} />
                {/each}
              {/each}

              <!-- Legend -->
              {#each d.lines as line, i}
                <circle cx={LINE_PAD.left + LINE_PW + 14} cy={LINE_PAD.top + 20 + i * 22} r="5"
                  fill={line.color} />
                <text x={LINE_PAD.left + LINE_PW + 22} y={LINE_PAD.top + 24 + i * 22}
                  font-size="10" fill="#475569" font-family="IBM Plex Mono, monospace">{line.habitat}</text>
              {/each}
            </svg>

          {:else}
            <div class="chart-placeholder">
              <div class="spinner"></div>
              <span>Loading chart…</span>
            </div>
          {/if}

        </div>
      </div>
    </div>
  {/if}

</div>

<style>
  .panel {
    display: flex;
    flex-direction: column;
    height: 100%;
    background: #f8fafc;
    overflow: hidden;
  }

  /* ── Study selector ──────────────────────────────────────────────────── */
  .study-selector {
    display: flex;
    gap: .5rem;
    padding: .85rem 1.25rem .75rem;
    border-bottom: 1px solid #e2e8f0;
    background: #fff;
    flex-shrink: 0;
  }

  .study-btn {
    padding: .38rem .85rem;
    border: 1px solid #e2e8f0;
    border-radius: 20px;
    background: #f8fafc;
    color: #64748b;
    font-size: .72rem;
    font-family: 'IBM Plex Sans', sans-serif;
    font-weight: 500;
    cursor: pointer;
    transition: all .15s;
    white-space: nowrap;
  }
  .study-btn:hover { border-color: #0854a8; color: #0854a8; }
  .study-btn.active {
    background: #0854a8;
    border-color: #0854a8;
    color: #fff;
  }

  /* ── Empty state ─────────────────────────────────────────────────────── */
  .empty {
    flex: 1;
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    gap: .75rem;
    color: #94a3b8;
    font-family: 'IBM Plex Sans', sans-serif;
    font-size: .85rem;
    text-align: center;
    padding: 2rem;
  }
  .empty p { max-width: 360px; line-height: 1.6; }
  .empty-icon { font-size: 2rem; color: #cbd5e1; }
  .err { color: #dc2626; font-family: 'IBM Plex Mono', monospace; font-size: .78rem; }

  .spinner {
    width: 28px; height: 28px;
    border: 2px solid #e2e8f0; border-top-color: #0854a8;
    border-radius: 50%; animation: spin .75s linear infinite;
  }
  @keyframes spin { to { transform: rotate(360deg); } }

  /* ── Study body ──────────────────────────────────────────────────────── */
  .study-body {
    flex: 1;
    overflow-y: auto;
    padding: 1.5rem 1.75rem;
  }

  .study-header {
    display: flex;
    flex-direction: column;
    gap: 1.5rem;
  }

  .study-meta {
    display: flex;
    flex-direction: column;
    gap: .5rem;
  }

  .study-title {
    font-size: 1.2rem;
    font-weight: 600;
    color: #0f172a;
    font-family: 'IBM Plex Sans', sans-serif;
    margin: 0;
  }

  .study-subtitle {
    font-size: .82rem;
    color: #0854a8;
    font-family: 'IBM Plex Mono', monospace;
    margin: 0;
  }

  .study-desc {
    font-size: .88rem;
    color: #475569;
    line-height: 1.7;
    margin: 0;
    font-family: 'IBM Plex Sans', sans-serif;
    max-width: 640px;
  }

  /* ── Chart ───────────────────────────────────────────────────────────── */
  .chart-area {
    background: #fff;
    border: 1px solid #e2e8f0;
    border-radius: 10px;
    padding: 1rem;
    box-shadow: 0 1px 4px rgba(0,0,0,.05);
  }

  .chart-svg {
    width: 100%;
    height: auto;
    display: block;
    font-family: 'IBM Plex Mono', monospace;
  }

  .chart-placeholder {
    display: flex;
    align-items: center;
    justify-content: center;
    gap: .75rem;
    height: 200px;
    color: #94a3b8;
    font-size: .8rem;
    font-family: 'IBM Plex Mono', monospace;
  }
</style>
