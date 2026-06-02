<script lang="ts">
  import { onMount } from 'svelte';
  import { db, runQuery, registerParquet } from '$lib/db.svelte';
  import * as Plot from '@observablehq/plot';

  type ExampleId = 'acanthurus' | 'uca' | 'gadus';
  type CellRow   = { cell: string; value: number; species?: string };

  let {
    activeExample,
    selectedSpecies,
    onSpeciesList,
    selectedProduct = 'glorysSST',
    selectedYear    = 'all',
    onYearList      = (_: string[]) => {},
  }: {
    activeExample:   'acanthurus' | 'uca' | 'gadus';
    selectedSpecies: string;
    onSpeciesList:   (list: string[]) => void;
    selectedProduct: string;
    selectedYear:    string;
    onYearList:      (list: string[]) => void;
  } = $props();

  let MapViewComp: any = $state(null);

  let mapRows     = $state<CellRow[]>([]);
  let mapLoading  = $state(false);
  let metricLabel = $state('');

  let chartData      = $state<any[]>([]);
  let chartLoading   = $state(false);
  let chartError     = $state('');
  let chartContainer = $state<HTMLDivElement | undefined>(undefined);

  let displayChartSQL = $state('');
  let displayMapSQL   = $state('');
  let activeSqlTab    = $state<'chart' | 'map'>('chart');
  let sqlCopied       = $state(false);

  let parquetLoaded = $state<Record<ExampleId, boolean>>({
    acanthurus: false, uca: false, gadus: false,
  });

  // ── Gadus constants ───────────────────────────────────────────────────────
  const GADUS_SST_DOMAIN  = ['glorysSST', 'coraltempSST', 'murSST', 'ostiaSST', 'glorysMaxDepthSST'];
  const GADUS_SST_COLORS  = ['#11b5ae',   '#4046ca',      '#f68512', '#de3c82',  '#7c3aed'];
  const GADUS_SST_LABELS: Record<string, string> = {
    glorysSST:         'GLORYS',
    coraltempSST:      'CoralTemp',
    murSST:            'MUR',
    ostiaSST:          'OSTIA',
    glorysMaxDepthSST: 'GLORYS\n(maximum\nrecord\ndepth)',
  };
  const GADUS_METRIC_LABELS: Record<string, string> = {
    glorysSST:         'avg GLORYS SST (°C)',
    coraltempSST:      'avg CoralTemp SST (°C)',
    murSST:            'avg MUR SST (°C)',
    ostiaSST:          'avg OSTIA SST (°C)',
    glorysMaxDepthSST: 'avg GLORYS depth temp (°C)',
  };

  // ── Canonical SQL for the SQL panel ──────────────────────────────────────
  function buildSQL(id: ExampleId) {
    if (id === 'acanthurus') {
      const cte = `WITH top_species AS (
  SELECT species
  FROM read_parquet('obistherm')
  WHERE family = 'Acanthuridae'
    AND absence IS NOT TRUE
    AND species IS NOT NULL
    AND deepTemperature IS NOT NULL
  GROUP BY species
  ORDER BY COUNT(*) DESC
  LIMIT 3
)`;
      return {
        chart: `${cte}
SELECT species, AphiaID, year, month, deepTemperature,
  h3_latlng_to_cell_string(decimalLatitude, decimalLongitude, 5) AS h3_5
FROM read_parquet('obistherm')
WHERE family = 'Acanthuridae'
  AND absence IS NOT TRUE
  AND species IN (SELECT species FROM top_species)
  AND deepTemperature IS NOT NULL`,
        map: `${cte}
SELECT h3_latlng_to_cell_string(decimalLatitude, decimalLongitude, 5) AS h3_5,
  AVG(deepTemperature) AS avg_deep_temp
FROM read_parquet('obistherm')
WHERE family = 'Acanthuridae'
  AND absence IS NOT TRUE
  AND species IN (SELECT species FROM top_species)
  AND deepTemperature IS NOT NULL
GROUP BY h3_5`,
      };
    }

    if (id === 'uca') {
      return {
        chart: `SELECT species, AphiaID, year, month,
  surfaceTemperature AS glorysSST,
  coraltempSST, murSST, ostiaSST,
  h3_latlng_to_cell_string(decimalLatitude, decimalLongitude, 5) AS h3_5
FROM read_parquet('obistherm')
WHERE species IN ('Minuca rapax', 'Leptuca thayeri')
  AND absence IS NOT TRUE`,
        map: `SELECT h3_latlng_to_cell_string(decimalLatitude, decimalLongitude, 5) AS h3_5,
  AVG(coraltempSST) AS avg_sst
FROM read_parquet('obistherm')
WHERE species IN ('Minuca rapax', 'Leptuca thayeri')
  AND absence IS NOT TRUE
  AND coraltempSST IS NOT NULL
GROUP BY h3_5`,
      };
    }

    // gadus
    return {
      chart: `SELECT
  h3_latlng_to_cell_string(decimalLatitude, decimalLongitude, 5) AS h3_5,
  MAKE_DATE(year, month, 1) AS date,
  AVG(surfaceTemperature)       AS glorysSST,
  AVG(coraltempSST)             AS coraltempSST,
  AVG(murSST)                   AS murSST,
  AVG(ostiaSST)                 AS ostiaSST,
  AVG(maximumDepthTemperature)  AS glorysMaxDepthSST,
  COUNT(*)                      AS n
FROM read_parquet('obistherm')
WHERE species = 'Gadus morhua'
  AND absence IS NOT TRUE
  AND year IS NOT NULL
  AND month IS NOT NULL
GROUP BY h3_5, date
ORDER BY date, h3_5`,
      map: `SELECT
  h3_latlng_to_cell_string(decimalLatitude, decimalLongitude, 5) AS h3_5,
  AVG(coraltempSST) AS avg_sst
FROM read_parquet('obistherm')
WHERE species = 'Gadus morhua'
  AND absence IS NOT TRUE
  AND coraltempSST IS NOT NULL
GROUP BY h3_5`,
    };
  }

  // ── Load data ─────────────────────────────────────────────────────────────
  async function loadData(
    id:      ExampleId,
    species: string,
    product: string,
    year:    string,
  ) {
    mapLoading   = true;
    chartLoading = true;
    chartError   = '';
    mapRows      = [];
    chartData    = [];

    const sql = buildSQL(id);
    displayChartSQL = sql.chart;
    displayMapSQL   = sql.map;

    try {
      const filename = `${id}.parquet`;

      if (!parquetLoaded[id]) {
        await registerParquet(filename, `/${filename}`);
        parquetLoaded[id] = true;

        const countRows = await runQuery(`SELECT COUNT(*) AS n FROM read_parquet('${filename}')`);
        db.rowCount += Number(countRows[0].n);
      }

      // Always repopulate species list — the page resets it on every example switch
      if (id === 'acanthurus') {
        const spRows = await runQuery(
          `SELECT DISTINCT species::VARCHAR AS species
           FROM read_parquet('acanthurus.parquet') ORDER BY species`
        );
        onSpeciesList(spRows.map(r => String(r.species)));
      } else if (id === 'uca') {
        onSpeciesList(['Minuca rapax', 'Leptuca thayeri']);
      }
      // gadus: year list is fetched below (must re-run whenever product changes)


      const spFilter = species !== 'all'
        ? `AND species = '${species.replace(/'/g, "''")}'`
        : '';

      // ── acanthurus ────────────────────────────────────────────────────────
      if (id === 'acanthurus') {
        const rows = await runQuery(`
          SELECT species::VARCHAR AS species,
                 year::INTEGER    AS year,
                 month::INTEGER   AS month,
                 AVG(deepTemperature)::DOUBLE AS avg_temp
          FROM read_parquet('acanthurus.parquet')
          WHERE deepTemperature IS NOT NULL ${spFilter}
          GROUP BY species, year, month
          ORDER BY species, year, month
        `);
        chartData = rows.map(r => ({
          species:  String(r.species),
          date:     new Date(Number(r.year), Number(r.month) - 1, 1),
          avg_temp: Number(r.avg_temp),
        }));

        const mRows = await runQuery(`
          SELECT h3_5::VARCHAR AS cell,
                 AVG(deepTemperature)::DOUBLE AS value,
                 STRING_AGG(DISTINCT species::VARCHAR, ', ') AS species
          FROM read_parquet('acanthurus.parquet')
          WHERE deepTemperature IS NOT NULL ${spFilter}
          GROUP BY h3_5
        `);
        mapRows     = mRows.map(r => ({ cell: String(r.cell), value: Number(r.value), species: r.species != null ? String(r.species) : undefined }));
        metricLabel = 'avg deepTemp (°C)';

      // ── uca ───────────────────────────────────────────────────────────────
      } else if (id === 'uca') {
        const rows = await runQuery(`
          SELECT species::VARCHAR     AS species,
                 year::INTEGER        AS year,
                 month::INTEGER       AS month,
                 glorysSST::DOUBLE    AS glorys,
                 coraltempSST::DOUBLE AS coraltemp,
                 murSST::DOUBLE       AS mur,
                 ostiaSST::DOUBLE     AS ostia
          FROM read_parquet('uca.parquet')
          WHERE 1=1 ${spFilter}
        `);

        const SOURCES = ['GLORYS', 'CoralTemp', 'MUR', 'OSTIA'] as const;
        const KEYS    = ['glorys', 'coraltemp', 'mur', 'ostia'] as const;
        const longData: any[] = [];
        for (const r of rows) {
          const date = new Date(Number(r.year), Number(r.month) - 1, 1);
          const sp   = String(r.species);
          for (let i = 0; i < SOURCES.length; i++) {
            const v = Number(r[KEYS[i]]);
            if (!isNaN(v)) longData.push({ species: sp, date, SSTsource: SOURCES[i], sst: v });
          }
        }
        chartData = longData;

        const mRows = await runQuery(`
          SELECT h3_5::VARCHAR AS cell,
                 AVG(coraltempSST)::DOUBLE AS value,
                 STRING_AGG(DISTINCT species::VARCHAR, ', ') AS species
          FROM read_parquet('uca.parquet')
          WHERE coraltempSST IS NOT NULL ${spFilter}
          GROUP BY h3_5
        `);
        mapRows     = mRows.map(r => ({ cell: String(r.cell), value: Number(r.value), species: r.species != null ? String(r.species) : undefined }));
        metricLabel = 'avg CoralTemp SST (°C)';

      // ── gadus ─────────────────────────────────────────────────────────────
      } else {
        // Fetch years that actually have data for the selected product.
        // This must run every time the product changes, not just on first load.
        const yRows = await runQuery(`
          SELECT DISTINCT EXTRACT(YEAR FROM date::DATE)::INTEGER AS year
          FROM read_parquet('gadus.parquet')
          WHERE ${product} IS NOT NULL
          ORDER BY year ASC
        `);
        const validYears = yRows.map(r => String(Number(r.year)));
        onYearList(validYears);

        // If the currently selected year has no data for this product, bail out
        // of the map query. The onYearList callback has already updated selectedYear
        // to the first valid year, which will trigger a clean re-run with correct data.
        const yearIsValid = year !== 'all' && validYears.includes(year);

        // Chart: yearly stats (avg/max/min) for all 5 SST products — always load
        const rows = await runQuery(`
          SELECT EXTRACT(YEAR FROM date::DATE)::INTEGER AS year,
            AVG(glorysSST)::DOUBLE         AS g_avg,
            MAX(glorysSST)::DOUBLE         AS g_max,
            MIN(glorysSST)::DOUBLE         AS g_min,
            AVG(coraltempSST)::DOUBLE      AS c_avg,
            MAX(coraltempSST)::DOUBLE      AS c_max,
            MIN(coraltempSST)::DOUBLE      AS c_min,
            AVG(murSST)::DOUBLE            AS m_avg,
            MAX(murSST)::DOUBLE            AS m_max,
            MIN(murSST)::DOUBLE            AS m_min,
            AVG(ostiaSST)::DOUBLE          AS o_avg,
            MAX(ostiaSST)::DOUBLE          AS o_max,
            MIN(ostiaSST)::DOUBLE          AS o_min,
            AVG(glorysMaxDepthSST)::DOUBLE AS d_avg,
            MAX(glorysMaxDepthSST)::DOUBLE AS d_max,
            MIN(glorysMaxDepthSST)::DOUBLE AS d_min
          FROM read_parquet('gadus.parquet')
          GROUP BY year
          ORDER BY year
        `);

        const SRC_MAP = [
          { key: 'glorysSST',         avg: 'g_avg', max: 'g_max', min: 'g_min' },
          { key: 'coraltempSST',       avg: 'c_avg', max: 'c_max', min: 'c_min' },
          { key: 'murSST',             avg: 'm_avg', max: 'm_max', min: 'm_min' },
          { key: 'ostiaSST',           avg: 'o_avg', max: 'o_max', min: 'o_min' },
          { key: 'glorysMaxDepthSST',  avg: 'd_avg', max: 'd_max', min: 'd_min' },
        ] as const;

        const longData: any[] = [];
        for (const r of rows) {
          const yr = Number(r.year);
          for (const s of SRC_MAP) {
            const rawAvg = r[s.avg];
            if (rawAvg !== null && rawAvg !== undefined) {
              longData.push({
                year:      yr,
                SSTsource: s.key,
                average:   Number(rawAvg),
                max:       Number(r[s.max]),
                min:       Number(r[s.min]),
              });
            }
          }
        }
        chartData = longData;

        // Map: only query when the year is confirmed valid for this product.
        // Use string prefix comparison on the date column to avoid DuckDB-WASM
        // predicate-pushdown issues with EXTRACT in WHERE on parquet date columns.
        if (yearIsValid) {
          const mapSql = `
            SELECT h3_5::VARCHAR AS cell,
                   AVG(${product})::DOUBLE AS value,
                   'Gadus morhua' AS species
            FROM read_parquet('gadus.parquet')
            WHERE ${product} IS NOT NULL
              AND CAST(date AS VARCHAR) LIKE '${year}%'
            GROUP BY h3_5
          `;
          const mRows = await runQuery(mapSql);
          mapRows     = mRows.map(r => ({ cell: String(r.cell), value: Number(r.value), species: String(r.species) }));
          metricLabel = GADUS_METRIC_LABELS[product] ?? `avg ${product} (°C)`;
        }
      }

    } catch (e) {
      chartError = String(e);
      console.error('ExamplesView loadData error:', e);
    } finally {
      mapLoading   = false;
      chartLoading = false;
    }
  }

  // ── Render chart ──────────────────────────────────────────────────────────
  $effect(() => {
    const data = chartData;
    const ex   = activeExample;
    const div  = chartContainer;
    if (!div) return;
    div.replaceChildren();
    if (!data.length) return;

    let chart: Element;

    if (ex === 'acanthurus') {
      chart = Plot.plot({
        width: 680, height: 460,
        marginLeft: 52, marginRight: 175, marginTop: 24, marginBottom: 36,
        marks: [
          Plot.frame({ fill: '#f8fafc', stroke: '#e2e8f0' }),
          Plot.line(data, { x: 'date', y: 'avg_temp', stroke: 'species', strokeOpacity: 0.25, fy: 'species' }),
          Plot.dot(data,  { x: 'date', y: 'avg_temp', fill: 'species', r: 2.5, fy: 'species' }),
        ],
        fy: { label: null },
        x:  { type: 'time', label: null, ticks: 5 },
        y:  { label: 'deepTemp (°C)', grid: true },
        color: { legend: true },
        style: { fontFamily: "'IBM Plex Mono', monospace", fontSize: '11px', background: 'transparent' },
      });

    } else if (ex === 'uca') {
      chart = Plot.plot({
        width: 620, height: 520,
        marginLeft: 52, marginRight: 85, marginTop: 24, marginBottom: 36,
        marks: [
          Plot.frame({ fill: '#f8fafc', stroke: '#e2e8f0' }),
          Plot.dot(data, { x: 'date', y: 'sst', fill: 'SSTsource', fillOpacity: 0.4, r: 2, fx: 'species', fy: 'SSTsource' }),
          Plot.linearRegressionY(data, { x: 'date', y: 'sst', stroke: '#334155', strokeWidth: 1.5, fx: 'species', fy: 'SSTsource' }),
        ],
        fx: { label: null },
        fy: { label: null, domain: ['GLORYS', 'CoralTemp', 'MUR', 'OSTIA'] },
        x:  { type: 'time', label: null, ticks: 4 },
        y:  { label: 'SST (°C)' },
        color: { domain: ['GLORYS', 'CoralTemp', 'MUR', 'OSTIA'], range: ['#11b5ae', '#4046ca', '#f68512', '#de3c82'], legend: false },
        style: { fontFamily: "'IBM Plex Mono', monospace", fontSize: '11px', background: 'transparent' },
      });

    } else {
      // gadus
      chart = Plot.plot({
        width: 620, height: 600,
        marginLeft: 52, marginRight: 115, marginTop: 24, marginBottom: 36,
        marks: [
          Plot.frame({ fill: '#f8fafc', stroke: '#e2e8f0' }),
          Plot.ruleX(data, {
            x: 'year', y1: 'min', y2: 'max',
            stroke: 'SSTsource', strokeWidth: 1.5, fy: 'SSTsource',
          }),
          Plot.dot(data, {
            x: 'year', y: 'average',
            fill: 'SSTsource', fillOpacity: 1, r: 3, fy: 'SSTsource',
          }),
        ],
        fy: {
          label: null,
          domain: GADUS_SST_DOMAIN,
          tickFormat: (d: string) => GADUS_SST_LABELS[d] ?? d,
        },
        x: { label: null, tickFormat: (d: number) => d.toString() },
        y: { label: 'Temperature (°C)' },
        color: { domain: GADUS_SST_DOMAIN, range: GADUS_SST_COLORS, legend: false },
        style: { fontFamily: "'IBM Plex Mono', monospace", fontSize: '11px', background: 'transparent' },
      });
    }

    div.appendChild(chart);
    return () => { div.replaceChildren(); };
  });

  // ── Reactive reload ───────────────────────────────────────────────────────
  $effect(() => {
    const id      = activeExample;
    const sp      = selectedSpecies;
    const product = selectedProduct;
    const year    = selectedYear;
    const status  = db.status;
    if (status !== 'ready') return;
    loadData(id, sp, product, year);
  });

  // ── SQL copy ──────────────────────────────────────────────────────────────
  function copySql() {
    navigator.clipboard.writeText(activeSqlTab === 'chart' ? displayChartSQL : displayMapSQL);
    sqlCopied = true;
    setTimeout(() => { sqlCopied = false; }, 1500);
  }

  onMount(async () => {
    MapViewComp = (await import('$lib/MapView.svelte')).default;
  });
</script>

<div class="wrap">

  <!-- ── Map panel ── -->
  <div class="map-panel">
    {#if MapViewComp}
      <MapViewComp rows={mapRows} loading={mapLoading} {metricLabel} bbox={null} />
    {:else}
      <div class="center-msg"><div class="spin"></div><span>Loading map…</span></div>
    {/if}
  </div>

  <!-- ── Chart panel ── -->
  <div class="chart-panel">

    <div class="chart-scroll">
      {#if chartLoading && !chartData.length}
        <div class="center-msg"><div class="spin"></div><span>Loading data…</span></div>
      {:else if chartError}
        <div class="center-msg err">{chartError}</div>
      {:else if !chartData.length}
        <div class="center-msg muted">Select an example to begin.</div>
      {:else}
        <p class="chart-title">
          {#if activeExample === 'acanthurus'}
            Acanthuridae — mean deep temperature over time (top 3 species)
          {:else if activeExample === 'uca'}
            Fiddler crabs — SST product comparison (<em>Minuca rapax</em> &amp; <em>Leptuca thayeri</em>)
          {:else}
            Atlantic cod (<em>Gadus morhua</em>) — temperature range by SST product and year
          {/if}
        </p>
        <div class="chart-inner" bind:this={chartContainer}></div>
      {/if}
    </div>

    {#if displayChartSQL}
      <div class="sql-panel">
        <div class="sql-header">
          <span class="sql-label">SQL</span>
          <div class="sql-tabs">
            <button class:active={activeSqlTab === 'chart'} onclick={() => { activeSqlTab = 'chart'; }}>Chart data</button>
            <button class:active={activeSqlTab === 'map'}   onclick={() => { activeSqlTab = 'map';   }}>Map data</button>
          </div>
          <button class="copy-btn" onclick={copySql}>{sqlCopied ? '✓ Copied' : 'Copy'}</button>
        </div>
        <textarea
          class="sql-box" readonly spellcheck="false"
          value={activeSqlTab === 'chart' ? displayChartSQL : displayMapSQL}
        ></textarea>
      </div>
    {/if}

  </div>
</div>

<style>
  .wrap { display: flex; height: 100%; overflow: hidden; }

  .map-panel {
    width: 44%; min-width: 280px;
    position: relative; overflow: hidden;
    border-right: 1px solid #e2e8f0;
  }

  .chart-panel {
    flex: 1; display: flex; flex-direction: column;
    overflow: hidden; background: #f8fafc;
  }

  .chart-scroll {
    flex: 1; overflow-y: auto;
    padding: 1.25rem 1.5rem 1rem;
    scrollbar-width: thin; scrollbar-color: rgba(0,0,0,.08) transparent;
  }

  .chart-title {
    font-size: .72rem; font-weight: 600; color: #334155;
    font-family: 'IBM Plex Mono', monospace;
    margin: 0 0 .85rem; letter-spacing: .01em;
  }

  .chart-inner :global(svg) { max-width: 100%; height: auto; display: block; }

  /* SQL panel */
  .sql-panel {
    height: 160px; min-height: 160px;
    border-top: 1px solid #e2e8f0; background: #ffffff;
    display: flex; flex-direction: column; flex-shrink: 0;
  }
  .sql-header {
    display: flex; align-items: center; gap: .6rem;
    padding: .4rem .85rem; border-bottom: 1px solid #e2e8f0; flex-shrink: 0;
  }
  .sql-label {
    font-family: 'IBM Plex Mono', monospace; font-size: .58rem; font-weight: 600;
    text-transform: uppercase; letter-spacing: .12em; color: #94a3b8; flex-shrink: 0;
  }
  .sql-tabs { display: flex; gap: .2rem; flex: 1; }
  .sql-tabs button {
    padding: .14rem .5rem; background: none; border: 1px solid transparent;
    border-radius: 4px; color: #94a3b8; font-size: .6rem;
    font-family: 'IBM Plex Mono', monospace; cursor: pointer; transition: all .15s;
  }
  .sql-tabs button:hover  { color: #475569; border-color: #e2e8f0; }
  .sql-tabs button.active { color: #0854a8; border-color: rgba(8,84,168,.3); background: rgba(8,84,168,.06); }
  .copy-btn {
    padding: .15rem .55rem; background: #fff; border: 1px solid #e2e8f0;
    border-radius: 4px; color: #64748b; font-size: .6rem;
    font-family: 'IBM Plex Mono', monospace; cursor: pointer; flex-shrink: 0;
  }
  .copy-btn:hover { color: #0854a8; border-color: #0854a8; }
  .sql-box {
    flex: 1; resize: none; background: transparent; border: none; outline: none;
    color: #0854a8; font-family: 'IBM Plex Mono', monospace;
    font-size: .72rem; line-height: 1.7; padding: .55rem .85rem; cursor: default;
    overflow-y: auto; scrollbar-width: thin; scrollbar-color: rgba(0,0,0,.08) transparent;
  }

  /* Helper */
  .center-msg {
    display: flex; align-items: center; justify-content: center;
    gap: .75rem; height: 200px; color: #94a3b8;
    font-size: .78rem; font-family: 'IBM Plex Mono', monospace;
  }
  .map-panel .center-msg { position: absolute; inset: 0; height: auto; }
  .center-msg.err   { color: #dc2626; }
  .center-msg.muted { color: #cbd5e1; }

  .spin {
    width: 20px; height: 20px; flex-shrink: 0;
    border: 2px solid #e2e8f0; border-top-color: #0854a8;
    border-radius: 50%; animation: spin .75s linear infinite;
  }
  @keyframes spin { to { transform: rotate(360deg); } }
</style>
