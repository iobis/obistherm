<script lang="ts">
  import { onMount } from 'svelte';
  import { browser } from '$app/environment';
  import { db, boot } from '$lib/db.svelte';
  import ExamplesView from '$lib/ExamplesView.svelte';
  import obislogo from '$lib/images/logo_simple.png';
  import ioclogo from '$lib/images/ioc_logo_black_2.svg';
  import { resolve } from '$app/paths';

  // ── Examples state ────────────────────────────────────────────────────────
  let activeExample      = $state<'acanthurus' | 'uca' | 'gadus'>('acanthurus');
  let selectedSpecies    = $state('all');
  let exampleSpeciesList = $state<string[]>([]);
  let selectedProduct    = $state('coraltempSST');
  let selectedYear       = $state('all');
  let gadusYearList      = $state<string[]>([]);

  // Reset all filters whenever the active example changes
  $effect(() => {
    void activeExample;
    selectedSpecies    = 'all';
    exampleSpeciesList = [];
    selectedProduct    = 'coraltempSST';
    selectedYear       = 'all';
    gadusYearList      = [];
  });

  function fmt(n: number) {
    if (n >= 1e6) return (n / 1e6).toFixed(1) + 'M';
    if (n >= 1e3) return (n / 1e3).toFixed(1) + 'k';
    return n.toLocaleString();
  }

  onMount(async () => {
    if (!browser) return;
    await boot();
  });
</script>

<svelte:head>
  <title>obistherm explorer</title>
  <link rel="stylesheet" href="https://unpkg.com/maplibre-gl@5.0.0/dist/maplibre-gl.css" />
  <link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:ital,wght@0,400;0,500;1,400&family=IBM+Plex+Sans:wght@400;500;600&display=swap" rel="stylesheet" />
</svelte:head>

<div class="layout">

  <!-- ════════════════ SIDEBAR ════════════════ -->
  <aside class="sidebar">

    <!-- Brand -->
    <div class="brand">
      <div class="brand-left">
        <span class="logo-text"><strong>obistherm</strong></span>
        <span class="logo-sub">explorer</span>
      </div>
      <a href={resolve('/docs/')} class="docs-link">Docs ↗</a>
    </div>

    <!-- Description -->
    <div class="description">
      OBIS occurrence records matched with sea surface temperature data from four different satellite products.
      <a href={resolve('/docs/')}>Want to know more?</a>
    </div>

    <!-- DB status -->
    <div class="db-status" data-s={db.status}>
      <span class="dot"></span>
      {#if db.status === 'booting'}  Connecting to DuckDB…
      {:else if db.status === 'ready'} <!--{fmt(db.rowCount)}-->Records loaded
      {:else if db.status === 'error'} {db.error}
      {:else} Idle
      {/if}
    </div>

    <!-- Section header -->
    <div class="section-hdr">Examples</div>

    <!-- Example list -->
    <div class="example-list">
      {#each ([
        { id: 'acanthurus', label: 'Acanthuridae',  sub: 'Deep temp · top 3 species'   },
        { id: 'uca',        label: 'Fiddler crabs', sub: 'SST products · 2 species'    },
        { id: 'gadus',      label: 'Atlantic cod',  sub: 'Gadus morhua · all products' },
      ] as const) as ex}
        <button
          class="ex-btn"
          class:active={activeExample === ex.id}
          onclick={() => { activeExample = ex.id; }}
        >
          <span class="ex-label">{ex.label}</span>
          <span class="ex-sub">{ex.sub}</span>
        </button>
      {/each}

      <!-- Species filter — acanthurus & uca -->
      {#if exampleSpeciesList.length > 0}
        <div class="species-filter">
          <label class="filter-label" for="sp-select">Species</label>
          <select id="sp-select" class="filter-select" bind:value={selectedSpecies}>
            <option value="all">All species</option>
            {#each exampleSpeciesList as sp}
              <option value={sp}>{sp}</option>
            {/each}
          </select>
        </div>
      {/if}

      <!-- Product + year filters — gadus -->
      {#if activeExample === 'gadus'}
        <div class="species-filter">
          <label class="filter-label" for="product-select">SST product</label>
          <select id="product-select" class="filter-select" bind:value={selectedProduct}>
            <option value="glorysSST">GLORYS</option>
            <option value="coraltempSST">CoralTemp</option>
            <option value="murSST">MUR</option>
            <option value="ostiaSST">OSTIA</option>
            <option value="glorysMaxDepthSST">GLORYS (depth)</option>
          </select>
        </div>
        <div class="species-filter">
          <label class="filter-label" for="year-select">Year</label>
          <select id="year-select" class="filter-select" bind:value={selectedYear}>
            <option value="all">All years</option>
            {#each gadusYearList as yr}
              <option value={yr}>{yr}</option>
            {/each}
          </select>
        </div>
      {/if}
    </div>

    <!-- Footer logos -->
    <div class="logo-footer">
      <img src={ioclogo} alt="IOC" style="height: 30px;" />
      <img src={obislogo} alt="OBIS" style="height: 27px;" />
    </div>

    <div class="sidebar-footer">
      <a href={resolve('/docs/')}>Documentation</a> ·
      <a href="https://obis.org" target="_blank" rel="noopener">obis.org</a> ·
      <a href="https://github.com/iobis/obis-therm" target="_blank" rel="noopener">GitHub</a>
    </div>

  </aside>

  <!-- ════════════════ MAIN ════════════════ -->
  <div class="main">
    <ExamplesView
      {activeExample}
      {selectedSpecies}
      {selectedProduct}
      {selectedYear}
      onSpeciesList={(list) => { exampleSpeciesList = list; }}
      onYearList={(list: string[]) => {
        gadusYearList = list;
        // Auto-select the first year that has data for this product.
        // Also re-selects if the user's current year has no data for a newly chosen product.
        if (list.length > 0 && (selectedYear === 'all' || !list.includes(selectedYear))) {
          selectedYear = list[0];
        }
      }}
    />
  </div>

</div>

<style>
  :global(*, *::before, *::after) { box-sizing: border-box; margin: 0; padding: 0; }
  :global(body) {
    background: #f0f4f8; color: #1e293b;
    font-family: 'IBM Plex Sans', sans-serif; overflow: hidden;
  }

  .layout { display: flex; height: 100vh; width: 100vw; }

  /* ── Sidebar ──────────────────────────────────────────────────────────── */
  .sidebar {
    width: 272px; min-width: 272px; height: 100vh;
    background: #ffffff;
    border-right: 1px solid #e2e8f0;
    display: flex; flex-direction: column;
    overflow-y: auto; overflow-x: hidden;
    scrollbar-width: thin; scrollbar-color: rgba(0,0,0,.08) transparent;
    font-family: 'IBM Plex Mono', monospace; font-size: .78rem;
  }

  .brand {
    display: flex; align-items: center; justify-content: space-between;
    padding: .85rem 1rem .7rem;
    border-bottom: 1px solid #e2e8f0;
  }
  .brand-left { display: flex; flex-direction: column; gap: .05rem; }
  .logo-text { font-size: .92rem; color: #0854a8; letter-spacing: .02em; }
  .logo-sub  { font-size: .58rem; color: #94a3b8; letter-spacing: .04em; text-transform: uppercase; }
  .docs-link { font-size: .62rem; color: #94a3b8; text-decoration: none; }
  .docs-link:hover { color: #0854a8; }

  .description {
    padding: .6rem 1rem .65rem;
    font-size: .67rem; color: #64748b; line-height: 1.55;
    font-family: 'IBM Plex Sans', sans-serif;
    border-bottom: 1px solid #f1f5f9;
  }
  .description a { color: #0854a8; text-decoration: none; white-space: nowrap; }
  .description a:hover { text-decoration: underline; }

  .db-status {
    display: flex; align-items: center; gap: .4rem;
    padding: .35rem 1rem; font-size: .62rem; color: #94a3b8;
    border-bottom: 1px solid #f1f5f9;
  }
  .db-status[data-s="ready"]   { color: #16a34a; }
  .db-status[data-s="booting"] { color: #0854a8; }
  .db-status[data-s="error"]   { color: #dc2626; }
  .dot { width: 5px; height: 5px; border-radius: 50%; background: currentColor; flex-shrink: 0; }
  .db-status[data-s="ready"] .dot { animation: blink 2.5s infinite; }
  @keyframes blink { 0%,100%{opacity:1} 50%{opacity:.2} }

  /* Section header */
  .section-hdr {
    padding: .48rem 1rem .4rem;
    font-size: .58rem; font-weight: 600; color: #94a3b8;
    text-transform: uppercase; letter-spacing: .1em;
    border-bottom: 1px solid #e2e8f0;
  }

  /* Example list */
  .example-list {
    display: flex; flex-direction: column;
    padding: .75rem; gap: .4rem; flex: 1;
  }
  .ex-btn {
    padding: .55rem .85rem; border: 1px solid #e2e8f0; border-radius: 7px;
    background: #f8fafc; color: #475569;
    font-size: .72rem; font-family: 'IBM Plex Mono', monospace;
    font-weight: 500; cursor: pointer; transition: all .15s; text-align: left;
    display: flex; flex-direction: column; gap: .18rem;
  }
  .ex-btn:hover  { border-color: #0854a8; color: #0854a8; background: rgba(8,84,168,.04); }
  .ex-btn.active { background: rgba(8,84,168,.08); border-color: rgba(8,84,168,.4); color: #0854a8; }
  .ex-label { font-weight: 600; font-size: .72rem; }
  .ex-sub   { font-size: .6rem; color: #94a3b8; font-weight: 400; }
  .ex-btn.active .ex-sub { color: rgba(8,84,168,.6); }

  /* Species filter */
  .species-filter {
    margin-top: .35rem; padding: .55rem .85rem .6rem;
    border: 1px solid #e2e8f0; border-radius: 7px; background: #fff;
    display: flex; flex-direction: column; gap: .3rem;
  }
  .filter-label {
    font-size: .58rem; font-weight: 600; color: #94a3b8;
    text-transform: uppercase; letter-spacing: .08em;
    font-family: 'IBM Plex Mono', monospace;
  }
  .filter-select {
    width: 100%; padding: .3rem .4rem;
    border: 1px solid #e2e8f0; border-radius: 5px;
    background: #f8fafc; color: #334155;
    font-size: .68rem; font-family: 'IBM Plex Mono', monospace;
    cursor: pointer; outline: none;
  }
  .filter-select:focus { border-color: #0854a8; }

  /* Footer */
  .logo-footer {
    display: flex; gap: .8rem;
    padding: .85rem 1rem;
    margin-top: auto;
  }
  .sidebar-footer {
    padding: .6rem 1rem .85rem;
    font-size: .58rem; color: #cbd5e1; line-height: 2;
    border-top: 1px solid #f1f5f9;
  }
  .sidebar-footer a { color: #94a3b8; text-decoration: none; }
  .sidebar-footer a:hover { color: #0854a8; }

  /* ── Main ────────────────────────────────────────────────────────────── */
  .main { flex: 1; overflow: hidden; }
</style>
