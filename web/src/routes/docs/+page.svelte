<script lang="ts">
  // Pure content page — no logic needed
</script>

<svelte:head>
  <title>obistherm — Documentation</title>
  <link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:ital,wght@0,400;0,500;0,600;1,400&display=swap" rel="stylesheet" />
</svelte:head>

<div class="page">

  <!-- Nav -->
  <nav>
    <a href="/" class="back">← Back to Explorer</a>
    <a href="https://github.com/iobis/obis-therm" target="_blank" rel="noopener" class="gh">
      GitHub ↗
    </a>
  </nav>

  <main>

    <!-- Hero -->
    <header class="hero">
      <div class="hero-eyebrow">IOC's Ocean Biodiversity Information System</div>
      <h1>obistherm</h1>
      <p class="hero-sub">
        OBIS occurrence records matched with multi-source sea temperature data.
      </p>
      <div class="hero-badges">
        <span class="badge blue">GeoParquet on S3</span>
        <span class="badge green">4 SST products</span>
        <span class="badge purple">1982–2025</span>
        <span class="badge orange">WoRMS aligned</span>
      </div>
    </header>

    <!-- What is it -->
    <section>
      <h2>What is obistherm?</h2>
      <p>
        <strong>obistherm</strong> is a dataset that extends species occurrence records from the 
        <a href="https://obis.org" target="_blank" rel="noopener">Ocean Biodiversity Information System (OBIS)</a> 
        with monthly seawater temperature data derived from four global satellite-based products. 
        Temperature values are matched to each occurrence location and, where available,
         to the recorded depth, with additional temperature estimates provided for the site's 
         maximum and mid-depth.
      </p>
      <p>
        The dataset is stored as <strong>GeoParquet</strong> on Amazon S3, making it queryable
        directly from the cloud using DuckDB, Python (arrow/pandas), or R — no download required.
        Each record retains the full taxonomic hierarchy aligned to WoRMS and the matched depth
        profile from GLORYS at surface, mid-water, and bottom levels.
      </p>

      <div class="callout blue">
        <div class="callout-icon">◈</div>
        <div>
          <strong>Serverless querying</strong> — This explorer uses DuckDB entirely in your
          browser (WebAssembly) using synthetic datasets. The real dataset on S3
          contains millions of records and can be queried with the same SQL shown here.
        </div>
      </div>
    </section>

    <!-- Temperature sources -->
    <section>
      <h2>Temperature sources</h2>
      <div class="card-grid">
        <div class="card">
          <div class="card-tag">GLORYS · CMEMS</div>
          <h3>GLORYS12</h3>
          <p>
            Global ocean reanalysis at 1/12° horizontal resolution (~8 km) with 50 vertical levels.
            Covers 1993 to present. Provides surface, mid-water, and bottom temperature profiles —
            the only product with subsurface data in obistherm.
          </p>
        </div>
        <div class="card">
          <div class="card-tag">CoralTemp · NOAA</div>
          <h3>CoralTemp Daily SST</h3>
          <p>
            Daily global 5 km satellite-derived SST from NOAA Coral Reef Watch.
            Coverage from 1985 to present. Designed for coral bleaching monitoring,
            widely used in marine ecology for its long record and global consistency.
          </p>
        </div>
        <div class="card">
          <div class="card-tag">MUR · NASA JPL</div>
          <h3>MUR SST (GHRSST L4)</h3>
          <p>
            Highest-resolution SST product in obistherm at 0.01° (~1 km), available from 2002.
            Multi-scale Ultra-high Resolution (MUR) blends multiple satellite and in-situ observations.
            Ideal for coastal and nearshore species analyses.
          </p>
        </div>
        <div class="card">
          <div class="card-tag">OSTIA · CMEMS/Met Office</div>
          <h3>OSTIA Foundation SST</h3>
          <p>
            Operational Sea Surface Temperature and Sea Ice Analysis at 0.05° resolution (~5 km),
            available from 2007. An L4 foundation SST product that removes the diurnal warming
            signal, making it suitable for longer-term climatological analyses.
          </p>
        </div>
      </div>
    </section>

    <!-- Dataset schema -->
    <section>
      <h2>Dataset schema</h2>
      <p>All columns in the obistherm dataset (Arrow schema types):</p>

      <div class="table-wrap">
        <table>
          <thead>
            <tr>
              <th>Column</th>
              <th>Type</th>
              <th>Description</th>
            </tr>
          </thead>
          <tbody>
            <tr><td><code>_id</code></td><td>string</td><td>Globally unique identifier assigned by OBIS</td></tr>
            <tr><td><code>dataset_id</code></td><td>string</td><td>Internal dataset identifier assigned by OBIS</td></tr>
            <tr><td><code>occurrenceID</code></td><td>string</td><td>Occurrence ID (DwC)</td></tr>
            <tr><td><code>datasetID</code></td><td>string</td><td>Dataset ID (DwC)</td></tr>
            <tr><td><code>AphiaID</code></td><td>int32</td><td>WoRMS AphiaID</td></tr>
            <tr><td><code>scientificName</code></td><td>string</td><td>Scientific name</td></tr>
            <tr><td><code>species</code></td><td>string</td><td>Species (from WoRMS)</td></tr>
            <tr><td><code>genus</code></td><td>string</td><td>Genus (from WoRMS)</td></tr>
            <tr><td><code>family</code></td><td>string</td><td>Family (from WoRMS)</td></tr>
            <tr><td><code>order</code></td><td>string</td><td>Order (from WoRMS)</td></tr>
            <tr><td><code>class</code></td><td>string</td><td>Class (from WoRMS)</td></tr>
            <tr><td><code>phylum</code></td><td>string</td><td>Phylum (from WoRMS)</td></tr>
            <tr><td><code>kingdom</code></td><td>string</td><td>Kingdom (from WoRMS)</td></tr>
            <tr><td><code>eventDate</code></td><td>string</td><td>Event date (DwC)</td></tr>
            <tr><td><code>date_start</code></td><td>double</td><td>Unix timestamp based on eventDate (start)</td></tr>
            <tr><td><code>date_mid</code></td><td>double</td><td>Unix timestamp based on eventDate (middle)</td></tr>
            <tr><td><code>date_end</code></td><td>double</td><td>Unix timestamp based on eventDate (end)</td></tr>
            <tr><td><code>decimalLongitude</code></td><td>double</td><td>Parsed and validated by OBIS</td></tr>
            <tr><td><code>decimalLatitude</code></td><td>double</td><td>Parsed and validated by OBIS</td></tr>
            <tr><td><code>coordinatePrecision</code></td><td>string</td><td>Precision of the coordinates</td></tr>
            <tr><td><code>coordinateUncertaintyInMeters</code></td><td>double</td><td>Uncertainty of the coordinates in meters</td></tr>
            <tr><td><code>minimumDepthInMeters</code></td><td>double</td><td>Maximum depth in meters</td></tr>
            <tr><td><code>maximumDepthInMeters</code></td><td>double</td><td>Minimum depth in meters</td></tr>
            <tr><td><code>absence</code></td><td>bool</td><td>If TRUE, is an absence</td></tr>
            <tr><td><code>flags</code></td><td>list&lt;string&gt;</td><td>OBIS QC flags</td></tr>
            <tr><td><code>month</code></td><td>int32</td><td>Month</td></tr>
            <tr><td><code>surfaceTemperature</code></td><td>double</td><td>GLORYS surface temperature</td></tr>
            <tr><td><code>midTemperature</code></td><td>double</td><td>GLORYS temperature at the mid-range of the water column</td></tr>
            <tr><td><code>deepTemperature</code></td><td>double</td><td>GLORYS temperature at the maximum range of depth</td></tr>
            <tr><td><code>bottomTemperature</code></td><td>double</td><td>GLORYS bottom temperature</td></tr>
            <tr><td><code>midDepth</code></td><td>double</td><td>Depth equivalent to the mid-range</td></tr>
            <tr><td><code>deepDepth</code></td><td>double</td><td>Depth equivalent to maximum depth</td></tr>
            <tr><td><code>minimumDepthTemperature</code></td><td>double</td><td>If <code>minimumDepthInMeters</code> is available, temperature at that depth</td></tr>
            <tr><td><code>maximumDepthTemperature</code></td><td>double</td><td>If <code>maximumDepthInMeters</code> is available, temperature at that depth</td></tr>
            <tr><td><code>minimumDepthClosestDepth</code></td><td>double</td><td>If the exact depth of <code>minimumDepthInMeters</code> is not available, the closest depth matched</td></tr>
            <tr><td><code>maximumDepthClosestDepth</code></td><td>double</td><td>If the exact depth of <code>maximumDepthInMeters</code> is not available, the closest depth matched</td></tr>
            <tr><td><code>coraltempSST</code></td><td>double</td><td>CoralTemp sea surface temperature</td></tr>
            <tr><td><code>murSST</code></td><td>double</td><td>MUR sea surface temperature</td></tr>
            <tr><td><code>ostiaSST</code></td><td>double</td><td>OSTIA sea surface temperature</td></tr>
            <tr><td><code>ostiaProduct</code></td><td>string</td><td>Which OSTIA product was used. REP = reprocessed, NRT = near real time</td></tr>
            <tr><td><code>obistherm_flags</code></td><td>string</td><td>obistherm flags</td></tr>
            <tr><td><code>h3_7</code></td><td>string</td><td>H3 grid cell at resolution 7</td></tr>
            <tr><td><code>geometry</code></td><td>binary</td><td>Geometry in WKB format</td></tr>
            <tr><td><code>year</code></td><td>int32</td><td>Year (partition column)</td></tr>
          </tbody>
        </table>
      </div>
    </section>

    <!-- Data access -->
    <section>
      <h2>Data access</h2>
      <p>The full dataset is available on AWS S3 as partitioned GeoParquet files:</p>

      <div class="code-block">
        <div class="code-label">S3 path</div>
        <pre><code>s3://obis-products/obistherm/</code></pre>
      </div>

      <div class="code-block">
        <div class="code-label">Download with AWS CLI (no credentials needed)</div>
        <pre><code>aws s3 cp --recursive s3://obis-products/obistherm . --no-sign-request</code></pre>
      </div>

      <div class="code-block">
        <div class="code-label">Query directly with DuckDB (Python)</div>
        <pre><code>import duckdb

con = duckdb.connect()
con.execute("INSTALL httpfs; LOAD httpfs; SET s3_region='us-east-1';")

# Mean SST by family
df = con.execute("""
    SELECT family,
           ROUND(AVG(coraltempSST), 2) AS mean_coraltemp,
           ROUND(AVG(murSST),       2) AS mean_mur,
           COUNT(*) AS records
    FROM read_parquet('s3://obis-products/obistherm/**/*.parquet')
    WHERE coraltempSST IS NOT NULL
    GROUP BY family
    ORDER BY mean_coraltemp DESC
""").df()
</code></pre>
      </div>

      <div class="code-block">
        <div class="code-label">Query with R + arrow</div>
        <pre><code>library(arrow)
library(dplyr)

ds &lt;- open_dataset("s3://obis-products/obistherm/",
                   filesystem = S3FileSystem$create(anonymous = TRUE))

ds |&gt;
  filter(!is.na(coraltempSST)) |&gt;
  group_by(family) |&gt;
  summarise(mean_sst = mean(coraltempSST), n = n()) |&gt;
  collect() |&gt;
  arrange(desc(mean_sst))
</code></pre>
      </div>

      <div class="code-block">
        <div class="code-label">Query with R + DuckDB</div>
        <pre><code>library(duckdb)
library(DBI)

con &lt;- dbConnect(duckdb())
dbExecute(con, "INSTALL httpfs; LOAD httpfs; SET s3_region='us-east-1';")

# Annual temperature trend
trend &lt;- dbGetQuery(con, "
  SELECT year,
         ROUND(AVG(coraltempSST), 3) AS coraltemp_mean,
         ROUND(AVG(murSST),       3) AS mur_mean,
         COUNT(*)                    AS records
  FROM read_parquet('s3://obis-products/obistherm/**/*.parquet')
  WHERE coraltempSST IS NOT NULL
  GROUP BY year
  ORDER BY year
")
</code></pre>
      </div>
    </section>

    <!-- Use cases -->
    <section>
      <h2>Use cases</h2>
      <div class="usecase-list">
        <div class="usecase">
          <h3>Thermal niche modeling</h3>
          <p>
            Match SST data directly to species occurrence records for species distribution
            modeling (SDM). The obistherm format eliminates the need to separately download
            and co-register SST rasters — thermal data is already paired with each observation.
          </p>
        </div>
        <div class="usecase">
          <h3>Climate change signals</h3>
          <p>
            Track how the temperatures experienced by marine species have changed over time.
            With records spanning 1982–2025 and four SST products, obistherm supports detection
            of long-term warming trends in species-level occurrence data.
          </p>
        </div>
        <div class="usecase">
          <h3>Multi-product comparison</h3>
          <p>
            Assess uncertainty from SST product choice in marine ecological analyses.
            Compare CoralTemp, MUR, OSTIA, and GLORYS at the same occurrence locations
            to understand product-level biases and their effect on thermal niche estimates.
          </p>
        </div>
      </div>
    </section>

    <!-- Citation -->
    <section>
      <h2>How to cite</h2>
      <p>If you use obistherm in your research, please cite:</p>
      <div class="citation-block">
        <p>Ocean Biodiversity Information System (OBIS) (25 March 2025) OBIS Occurrence Data. Ocean Biodiversity Information System. Intergovernmental Oceanographic Commission of UNESCO. https://obis.org.</p>
        <p>Copernicus Marine Environment Monitoring Service (CMEMS). GLORYS12V1 Global Ocean Physics Reanalysis. https://doi.org/10.48670/moi-00021</p>
        <p>NOAA Coral Reef Watch. CoralTemp Daily Global 5km Satellite Coral Bleaching Monitoring Products. https://coralreefwatch.noaa.gov</p>
        <p>JPL MUR MEaSUREs Project (2015). GHRSST Level 4 MUR Global Foundation Sea Surface Temperature Analysis (v4.1). NASA Physical Oceanography DAAC. https://doi.org/10.5067/GHGMR-4FJ04</p>
      </div>
    </section>

  </main>
</div>

<style>
  :global(*, *::before, *::after) { box-sizing: border-box; margin: 0; padding: 0; }
  :global(body) {
    background: #f0f4f8;
    color: #334155;
    font-family: 'IBM Plex Sans', sans-serif;
    overflow-y: auto;
  }

  .page { min-height: 100vh; }

  nav {
    position: sticky; top: 0; z-index: 10;
    display: flex; align-items: center; justify-content: space-between;
    padding: 0.75rem 2rem;
    background: rgba(255,255,255,0.95);
    border-bottom: 1px solid #e2e8f0;
    backdrop-filter: blur(8px);
  }
  .back { font-size: 0.78rem; color: #64748b; text-decoration: none; transition: color 0.15s; }
  .back:hover { color: #0854a8; }
  .gh   { font-size: 0.72rem; color: #94a3b8; text-decoration: none; transition: color 0.15s; }
  .gh:hover { color: #0854a8; }

  main {
    max-width: 800px;
    margin: 0 auto;
    padding: 3rem 2rem 5rem;
    display: flex;
    flex-direction: column;
    gap: 3.5rem;
  }

  /* ── Hero ────────────────────────────────────────────────────────────── */
  .hero { display: flex; flex-direction: column; gap: 1rem; }
  .hero-eyebrow {
    font-size: 0.7rem; font-weight: 500;
    text-transform: uppercase; letter-spacing: 0.14em;
    color: #0854a8;
  }
  h1 {
    font-size: 3rem; font-weight: 600; color: #0f172a;
    line-height: 1.1; letter-spacing: -0.02em;
    font-family: 'IBM Plex Sans', sans-serif;
  }
  .hero-sub {
    font-size: 1.05rem; color: #64748b;
    line-height: 1.6; max-width: 560px;
  }
  .hero-badges { display: flex; flex-wrap: wrap; gap: 0.5rem; margin-top: 0.5rem; }
  .badge {
    font-family: 'IBM Plex Mono', monospace;
    font-size: 0.65rem; font-weight: 500;
    padding: 0.2rem 0.65rem;
    border-radius: 20px; border: 1px solid;
  }
  .badge.blue   { color: #0854a8; border-color: rgba(8,84,168,0.25);  background: rgba(8,84,168,0.06);  }
  .badge.green  { color: #16a34a; border-color: rgba(22,163,74,0.25); background: rgba(22,163,74,0.06); }
  .badge.purple { color: #7c3aed; border-color: rgba(124,58,237,0.25);background: rgba(124,58,237,0.06);}
  .badge.orange { color: #ea580c; border-color: rgba(234,88,12,0.25); background: rgba(234,88,12,0.06); }

  /* ── Sections ────────────────────────────────────────────────────────── */
  section { display: flex; flex-direction: column; gap: 1.25rem; }
  h2 {
    font-size: 1.35rem; font-weight: 600; color: #0f172a;
    border-bottom: 1px solid #e2e8f0;
    padding-bottom: 0.6rem;
  }
  p { font-size: 0.92rem; color: #475569; line-height: 1.75; }
  a { color: #0854a8; text-decoration: none; }
  a:hover { text-decoration: underline; }

  /* ── Callout ─────────────────────────────────────────────────────────── */
  .callout {
    display: flex; gap: 1rem; align-items: flex-start;
    padding: 1rem 1.25rem; border-radius: 8px;
    font-size: 0.88rem; color: #475569; line-height: 1.65;
  }
  .callout.blue {
    background: rgba(8,84,168,0.05);
    border: 1px solid rgba(8,84,168,0.18);
  }
  .callout-icon { font-size: 1.2rem; color: #0854a8; flex-shrink: 0; margin-top: 0.1rem; }
  .callout strong { color: #0f172a; }

  /* ── Cards ───────────────────────────────────────────────────────────── */
  .card-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 1rem; }
  .card {
    background: #fff; border: 1px solid #e2e8f0; border-radius: 8px;
    padding: 1.1rem; display: flex; flex-direction: column; gap: 0.5rem;
    box-shadow: 0 1px 3px rgba(0,0,0,0.05);
  }
  .card-tag {
    font-family: 'IBM Plex Mono', monospace;
    font-size: 0.62rem; font-weight: 500;
    color: #0854a8; letter-spacing: 0.05em; text-transform: uppercase;
  }
  .card h3 { font-size: 0.9rem; font-weight: 600; color: #0f172a; }
  .card p   { font-size: 0.8rem; color: #64748b; line-height: 1.65; }

  /* ── Table ───────────────────────────────────────────────────────────── */
  .table-wrap {
    overflow-x: auto; border-radius: 8px;
    border: 1px solid #e2e8f0; background: #fff;
    box-shadow: 0 1px 3px rgba(0,0,0,0.05);
  }
  table { width: 100%; border-collapse: collapse; font-size: 0.82rem; }
  th {
    background: #f8fafc; color: #64748b; font-weight: 600; font-size: 0.72rem;
    text-transform: uppercase; letter-spacing: 0.07em;
    padding: 0.6rem 0.85rem; text-align: left; border-bottom: 1px solid #e2e8f0;
  }
  td { padding: 0.5rem 0.85rem; color: #475569; border-bottom: 1px solid #f1f5f9; }
  tr:last-child td { border-bottom: none; }
  tr:hover td { background: #f8fafc; }
  code {
    font-family: 'IBM Plex Mono', monospace; font-size: 0.78em;
    color: #0854a8; background: rgba(8,84,168,0.07);
    padding: 0.1rem 0.35rem; border-radius: 3px;
  }

  /* ── Code blocks ─────────────────────────────────────────────────────── */
  .code-block {
    border-radius: 8px; overflow: hidden;
    border: 1px solid #e2e8f0; box-shadow: 0 1px 3px rgba(0,0,0,0.05);
  }
  .code-label {
    background: #f8fafc; padding: 0.35rem 0.85rem;
    font-size: 0.62rem; color: #64748b;
    font-family: 'IBM Plex Mono', monospace;
    letter-spacing: 0.05em; border-bottom: 1px solid #e2e8f0;
  }
  pre { background: #0f172a; padding: 0.85rem 1rem; overflow-x: auto; margin: 0; }
  pre code {
    font-family: 'IBM Plex Mono', monospace; font-size: 0.78rem;
    color: #7dd3fc; background: none; padding: 0; border-radius: 0;
    line-height: 1.65; display: block;
  }

  /* ── Use cases ───────────────────────────────────────────────────────── */
  .usecase-list { display: flex; flex-direction: column; gap: 1.25rem; }
  .usecase {
    padding: 1rem 1.25rem;
    border-left: 3px solid rgba(8,84,168,0.3);
    background: #fff; border-radius: 0 7px 7px 0;
    box-shadow: 0 1px 3px rgba(0,0,0,0.04);
  }
  .usecase h3 { font-size: 0.9rem; font-weight: 600; color: #0f172a; margin-bottom: 0.4rem; }
  .usecase p  { font-size: 0.82rem; color: #64748b; line-height: 1.65; }

  /* ── Citation ────────────────────────────────────────────────────────── */
  .citation-block {
    background: #fff; border: 1px solid #e2e8f0; border-radius: 7px;
    padding: 1rem 1.25rem; display: flex; flex-direction: column; gap: 0.5rem;
    box-shadow: 0 1px 3px rgba(0,0,0,0.04);
  }
  .citation-block p {
    font-family: 'IBM Plex Mono', monospace;
    font-size: 0.72rem; color: #64748b; line-height: 1.6;
  }
</style>
