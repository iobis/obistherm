import * as duckdb from '@duckdb/duckdb-wasm';
import { latLngToCell } from 'h3-js';

export const db = $state({
  status: 'idle' as 'idle' | 'booting' | 'ready' | 'error',
  error: null as string | null,
  rowCount: 0,
});

export let connection: duckdb.AsyncDuckDBConnection | null = null;
let duckInstance: duckdb.AsyncDuckDB | null = null;
let bootPromise: Promise<void> | null = null;

export async function boot(): Promise<void> {
  if (bootPromise) return bootPromise;
  bootPromise = _boot();
  return bootPromise;
}


async function _boot() {
  db.status = 'booting';
  try {
    const bundle = await duckdb.selectBundle(duckdb.getJsDelivrBundles());
    const workerUrl = URL.createObjectURL(
      new Blob([`importScripts("${bundle.mainWorker}");`], { type: 'text/javascript' })
    );
    const worker = new Worker(workerUrl);
    URL.revokeObjectURL(workerUrl);
    const instance = new duckdb.AsyncDuckDB(
      new duckdb.ConsoleLogger(duckdb.LogLevel.WARNING), worker
    );
    await instance.instantiate(bundle.mainModule, bundle.pthreadWorker);
    duckInstance = instance;
    connection = await instance.connect();

    await connection.query('INSTALL httpfs; LOAD httpfs;');

    db.status = 'ready';
  } catch (err) {
    db.error  = String(err);
    db.status = 'error';
    throw err;
  }
}

export async function runQuery(sql: string): Promise<Record<string, unknown>[]> {
  if (!connection) throw new Error('DuckDB not ready');
  const result = await connection.query(sql);
  const cols   = result.schema.fields.map((f: { name: string }) => f.name);
  return result.toArray().map((row: Record<string, unknown>) =>
    Object.fromEntries(cols.map((c: string) => [c, row[c]]))
  );
}

// Fetches a parquet file from `path` and registers it with DuckDB under `name`.
// After calling this, the file is queryable as: read_parquet('name')
export async function registerParquet(name: string, path: string): Promise<void> {
  if (!duckInstance) throw new Error('DuckDB not initialized');
  const res = await fetch(path);
  if (!res.ok) throw new Error(`Failed to fetch parquet ${path}: ${res.status}`);
  await duckInstance.registerFileBuffer(name, new Uint8Array(await res.arrayBuffer()));
}

// Returns {cell, value} for map display. SQL must return decimalLatitude, decimalLongitude, value
export async function runMapQuery(sql: string, h3Res = 5): Promise<{cell: string; value: number}[]> {
  const rows = await runQuery(sql);
  const cellMap = new Map<string, {sum: number; count: number}>();
  for (const row of rows) {
    const lat = Number(row.decimalLatitude);
    const lng = Number(row.decimalLongitude);
    const val = Number(row.value);
    if (isNaN(lat) || isNaN(lng) || isNaN(val)) continue;
    const cell = latLngToCell(lat, lng, h3Res);
    const e = cellMap.get(cell) ?? { sum: 0, count: 0 };
    e.sum += val; e.count++;
    cellMap.set(cell, e);
  }
  return Array.from(cellMap.entries()).map(([cell, { sum, count }]) => ({
    cell,
    value: sum / count,
  }));
}
