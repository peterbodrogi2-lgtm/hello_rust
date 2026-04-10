#![allow(non_snake_case)]
// ─── Crates ───────────────────────────────────────────────────────────────────
use axum::{
    extract::{Query, State},
    routing::get,
    Json, Router,
};
use csv::ReaderBuilder;
use serde::{Deserialize, Serialize};
use std::sync::Arc;
use tokio::net::TcpListener;

// ─── 2-D Look-Up Table ───────────────────────────────────────────────────────
//
//  Regular grid: I  = 0, 1, …, 2000 mA  (N_I  = 2001 rows)
//                Ts = 25, 26, …, 85 °C  (N_TS =   61 cols)
//  Storage: row-major flat Vec<f64>, index = i_idx * N_TS + ts_idx
//
struct Lut2D {
    n_i:    usize,
    n_ts:   usize,
    i_min:  f64,
    i_step: f64,   // = 1.0
    ts_min: f64,
    ts_step: f64,  // = 1.0
    data:   Vec<f64>,
}

impl Lut2D {
    /// Load a headerless comma-separated CSV.
    /// Each row is one I-value; each column is one Ts-value.
    fn load(path: &str, n_i: usize, n_ts: usize,
        i_min: f64, i_step: f64,
        ts_min: f64, ts_step: f64) -> anyhow::Result<Self> {
    let mut rdr = ReaderBuilder::new()
        .has_headers(true)   // ← consumes the 61-field Ts header row, discards it
        .flexible(true)      // ← allows data rows (62 fields) to differ from header (61 fields)
        .from_path(path)?;

    let mut data = Vec::with_capacity(n_i * n_ts);
    for rec in rdr.records() {
        let rec = rec?;
        // rec[0] is the R row-name (I value as string: "0","1",…,"2000") → skip it
        // rec[1..=61] are the 61 flux/power values for Ts = 25…85
        for field in rec.iter().skip(1).take(n_ts) {
            data.push(field.trim().parse::<f64>()
                .map_err(|e| anyhow::anyhow!("parse error in {path}: {e}"))?);
        }
    }
    anyhow::ensure!(
        data.len() == n_i * n_ts,
        "LUT '{path}': expected {} values, got {}",
        n_i * n_ts, data.len()
    );
    Ok(Lut2D { n_i, n_ts, i_min, i_step, ts_min, ts_step, data })
}

    /// Bilinear interpolation; inputs are clamped to the grid boundary.
    /// Equivalent to R's itp2df().
    #[inline]
    fn interp(&self, i: f64, ts: f64) -> f64 {
        let i_max  = self.i_min  + (self.n_i  - 1) as f64 * self.i_step;
        let ts_max = self.ts_min + (self.n_ts - 1) as f64 * self.ts_step;
        let i  = i .clamp(self.i_min,  i_max);
        let ts = ts.clamp(self.ts_min, ts_max);

        // Lower index in each axis, clamped so i1+1 is always valid
        let i1 = (((i  - self.i_min)  / self.i_step)  .floor() as usize).min(self.n_i  - 2);
        let t1 = (((ts - self.ts_min) / self.ts_step) .floor() as usize).min(self.n_ts - 2);
        let i2 = i1 + 1;
        let t2 = t1 + 1;

        // Fractional distances in [0,1]
        let di = (i  - (self.i_min  + i1 as f64 * self.i_step))  / self.i_step;
        let dt = (ts - (self.ts_min + t1 as f64 * self.ts_step)) / self.ts_step;

        let f11 = self.data[i1 * self.n_ts + t1];
        let f21 = self.data[i2 * self.n_ts + t1];
        let f12 = self.data[i1 * self.n_ts + t2];
        let f22 = self.data[i2 * self.n_ts + t2];

        f11 * (1.0 - di) * (1.0 - dt)
            + f21 * di       * (1.0 - dt)
            + f12 * (1.0 - di) * dt
            + f22 * di       * dt
    }
}

// ─── Brent's bounded scalar minimiser ────────────────────────────────────────
//
//  Direct Rust port of R's optimize() which uses Brent (1973).
//  Finds x* in [a, b] that minimises f(x).
//  Tolerance matches R's default (tol = .Machine$double.eps^0.25 ≈ 1.2e-4),
//  tightened here to 1e-8 for accuracy.
//
fn brent_min<F: Fn(f64) -> f64>(f: F, mut a: f64, mut b: f64) -> f64 {
    const C: f64 = 0.381_966_011_250_105; // (3 - √5) / 2  — golden ratio conjugate
    const TOL: f64 = 1e-8;

    let mut x = a + C * (b - a);
    let (mut w, mut v) = (x, x);
    let mut fx = f(x);
    let (mut fw, mut fv) = (fx, fx);
    let (mut d, mut e) = (0.0_f64, 0.0_f64);

    for _ in 0..500 {
        let m    = 0.5 * (a + b);
        let tol1 = TOL * x.abs() + 1e-10;
        let tol2 = 2.0 * tol1;

        // Convergence test
        if (x - m).abs() <= tol2 - 0.5 * (b - a) {
            return x;
        }

        // Try parabolic interpolation
        if e.abs() > tol1 {
            let r  = (x - w) * (fx - fv);
            let q  = (x - v) * (fx - fw);
            let mut p  = (x - v) * q - (x - w) * r;
            let mut q2 = 2.0 * (q - r);
            if q2 > 0.0 { p = -p; } else { q2 = -q2; }

            if p.abs() < (0.5 * q2 * e).abs()
                && p > q2 * (a - x)
                && p < q2 * (b - x)
            {
                // Parabolic step
                e = d;
                d = p / q2;
                let u_try = x + d;
                if (u_try - a) < tol2 || (b - u_try) < tol2 {
                    d = if x < m { tol1 } else { -tol1 };
                }
            } else {
                // Golden-section step
                e = if x < m { b - x } else { a - x };
                d = C * e;
            }
        } else {
            // Golden-section step
            e = if x < m { b - x } else { a - x };
            d = C * e;
        }

        let u  = x + if d.abs() >= tol1 { d } else if d >= 0.0 { tol1 } else { -tol1 };
        let fu = f(u);

        if fu <= fx {
            if u < x { b = x; } else { a = x; }
            v = w; fv = fw;
            w = x; fw = fx;
            x = u; fx = fu;
        } else {
            if u < x { a = u; } else { b = u; }
            if fu <= fw || w == x       { v = w; fv = fw; w = u; fw = fu; }
            else if fu <= fv || v == x || v == w { v = u; fv = fu; }
        }
    }
    x // return best found if max iterations reached
}

// ─── LED catalogue & LUT dimensions ─────────────────────────────────────────
const LED_NAMES: [&str; 5] = [
    "O_HP_2700_92",
    "O_HP_3000_92",
    "O_HP_3500_92",
    "O_HP_4000_92",
    "O_HP_3000_97",
];
const N_I:    usize = 2001; // rows: I = 0, 1, …, 2000 mA
const N_TS:   usize =   61; // cols: Ts = 25, 26, …, 85 °C
const I_MIN:  f64   =  0.0;
const I_STEP: f64   =  1.0;
const TS_MIN: f64   = 25.0;
const TS_STEP:f64   =  1.0;

// ─── Shared application state (LUTs loaded once at startup) ──────────────────
struct AppState {
    flux_luts:    [Lut2D; 5],
    thermal_luts: [Lut2D; 5],
}

// AppState is read-only after construction → Send + Sync via Arc
unsafe impl Send for AppState {}
unsafe impl Sync for AppState {}

// ─── HTTP query parameter struct ─────────────────────────────────────────────
//
//  Field names match the R/plumber API exactly (kept non-snake-case).
//  All parameters have defaults identical to the plumber.R defaults.
//
#[derive(Deserialize)]
struct ThermalParams {
    #[serde(default = "d_n_led")]      n_LED:        f64,
    #[serde(default = "d_tq_from")]    Tq_from:      f64,
    #[serde(default = "d_tq_to")]      Tq_to:        f64,
    #[serde(default = "d_tq_by")]      Tq_by:        f64,
    #[serde(default = "d_tq_ref")]     Tq_ref:       f64,
    #[serde(default = "d_tpcb")]       TPCBgrenz:    f64,
    #[serde(default = "d_iled")]       ILED_abs_max: f64,
    #[serde(default = "d_rpcb")]       RPCB_U:       f64,
    #[serde(default = "d_tu_ist")]     TU_ist:       f64,
}
fn d_n_led()   -> f64 { 16.0 }
fn d_tq_from() -> f64 { 20.0 }
fn d_tq_to()   -> f64 { 90.0 }
fn d_tq_by()   -> f64 {  5.0 }
fn d_tq_ref()  -> f64 { 25.0 }
fn d_tpcb()    -> f64 { 87.534_1 }
fn d_iled()    -> f64 { 870.0 }
fn d_rpcb()    -> f64 { 1.528_882_045_951 }
fn d_tu_ist()  -> f64 { 22.0 }

// ─── JSON response row ────────────────────────────────────────────────────────
#[derive(Serialize, Clone)]
struct ThermalRow {
    LED_name:            String,
    Tq:                  f64,
    max_allowed_current: f64,
    TPCBist:             f64,
    flux:                f64,
    thermal_power:       f64,
    TNTC_ist:            f64,
    rel_flux_in_percent: f64,
}

// ─── Core computation — equivalent to thermal_regulation_internal() in R ─────
fn regulate_one(
    fl:           &Lut2D,
    th:           &Lut2D,
    n_LED:        f64,
    Tq:           f64,
    TPCBgrenz:    f64,
    ILED_abs_max: f64,
    RPCB_U:       f64,
    TU_ist:       f64,
    led_name:     &str,
) -> ThermalRow {
    // ① Find NTC temperature at ILED_abs_max  (fixed-point solve)
    //    TNTC = TU_ist + RPCB_U * n_LED * P_thermal(ILED_abs_max, TNTC)
    let TNTC_ist = brent_min(
        |x| (x - (TU_ist + RPCB_U * n_LED * th.interp(ILED_abs_max, x))).abs(),
        0.0, 100.0,
    );

    // ② Allowed thermal power budget at Tq
    let PW: f64 = (TPCBgrenz - Tq) / RPCB_U;
    // TPCBist = Tq + PW * RPCB_U = TPCBgrenz (algebraically exact)
    let mut TPCBist: f64 = TPCBgrenz;

    // ③ Maximum current that produces exactly PW at TPCBist
    let mut max_allowed_current = brent_min(
        |x| (n_LED * th.interp(x, TPCBist) - PW).abs(),
        0.0, 2000.0,
    );

    // ④ If current exceeds absolute max: clamp and re-solve TPCBist
    if max_allowed_current > ILED_abs_max {
        max_allowed_current = ILED_abs_max;
        TPCBist = brent_min(
            |x| (x - (Tq + n_LED * th.interp(ILED_abs_max, x) * RPCB_U)).abs(),
            0.0, 100.0,
        );
    }

    let flux          = n_LED * fl.interp(max_allowed_current, TPCBist);
    let thermal_power = n_LED * th.interp(max_allowed_current, TPCBist);

    ThermalRow {
        LED_name: led_name.to_owned(),
        Tq,
        max_allowed_current,
        TPCBist,
        flux,
        thermal_power,
        TNTC_ist,
        rel_flux_in_percent: 0.0, // filled by caller
    }
}

// ─── Main regulation function — equivalent to thermal_regulation() in R ───────
fn thermal_regulation(state: &AppState, p: &ThermalParams) -> Vec<ThermalRow> {
    // Reference flux per LED at Tq_ref
    let ref_flux: Vec<f64> = (0..LED_NAMES.len())
        .map(|i| regulate_one(
            &state.flux_luts[i], &state.thermal_luts[i],
            p.n_LED, p.Tq_ref, p.TPCBgrenz, p.ILED_abs_max, p.RPCB_U, p.TU_ist,
            LED_NAMES[i],
        ).flux)
        .collect();

    // Build Tq sequence
    let n_steps = ((p.Tq_to - p.Tq_from) / p.Tq_by).round() as usize + 1;
    let mut rows: Vec<ThermalRow> = Vec::with_capacity(n_steps * LED_NAMES.len());

    for step in 0..n_steps {
        let tq = p.Tq_from + step as f64 * p.Tq_by;
        for (i, &name) in LED_NAMES.iter().enumerate() {
            let mut row = regulate_one(
                &state.flux_luts[i], &state.thermal_luts[i],
                p.n_LED, tq, p.TPCBgrenz, p.ILED_abs_max, p.RPCB_U, p.TU_ist,
                name,
            );
            row.rel_flux_in_percent = 100.0 * row.flux / ref_flux[i];
            rows.push(row);
        }
    }

    // Sort by LED_name (matches R: df_out3[order(df_out3$LED_name), ])
    rows.sort_by(|a, b| a.LED_name.cmp(&b.LED_name)
                         .then(a.Tq.partial_cmp(&b.Tq).unwrap()));
    rows
}

// ─── Axum GET handler ─────────────────────────────────────────────────────────
async fn echo_handler(
    State(state): State<Arc<AppState>>,
    Query(params): Query<ThermalParams>,
) -> Json<Vec<ThermalRow>> {
    // Offload CPU-bound work to Tokio's blocking thread pool
    // so the async executor stays responsive under concurrent requests.
    let result = tokio::task::spawn_blocking(move || {
        thermal_regulation(&state, &params)
    })
    .await
    .expect("thermal_regulation panicked");

    Json(result)
}

// ─── Entry point ─────────────────────────────────────────────────────────────
#[tokio::main]
async fn main() -> anyhow::Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| "info".into()),
        )
        .init();

    tracing::info!("Loading LUTs …");
    let flux_luts = [
        Lut2D::load(&format!("luts/flux_{}.csv",          LED_NAMES[0]), N_I, N_TS, I_MIN, I_STEP, TS_MIN, TS_STEP)?,
        Lut2D::load(&format!("luts/flux_{}.csv",          LED_NAMES[1]), N_I, N_TS, I_MIN, I_STEP, TS_MIN, TS_STEP)?,
        Lut2D::load(&format!("luts/flux_{}.csv",          LED_NAMES[2]), N_I, N_TS, I_MIN, I_STEP, TS_MIN, TS_STEP)?,
        Lut2D::load(&format!("luts/flux_{}.csv",          LED_NAMES[3]), N_I, N_TS, I_MIN, I_STEP, TS_MIN, TS_STEP)?,
        Lut2D::load(&format!("luts/flux_{}.csv",          LED_NAMES[4]), N_I, N_TS, I_MIN, I_STEP, TS_MIN, TS_STEP)?,
    ];
    let thermal_luts = [
        Lut2D::load(&format!("luts/thermal_power_{}.csv", LED_NAMES[0]), N_I, N_TS, I_MIN, I_STEP, TS_MIN, TS_STEP)?,
        Lut2D::load(&format!("luts/thermal_power_{}.csv", LED_NAMES[1]), N_I, N_TS, I_MIN, I_STEP, TS_MIN, TS_STEP)?,
        Lut2D::load(&format!("luts/thermal_power_{}.csv", LED_NAMES[2]), N_I, N_TS, I_MIN, I_STEP, TS_MIN, TS_STEP)?,
        Lut2D::load(&format!("luts/thermal_power_{}.csv", LED_NAMES[3]), N_I, N_TS, I_MIN, I_STEP, TS_MIN, TS_STEP)?,
        Lut2D::load(&format!("luts/thermal_power_{}.csv", LED_NAMES[4]), N_I, N_TS, I_MIN, I_STEP, TS_MIN, TS_STEP)?,
    ];
    tracing::info!("LUTs loaded OK — {} entries each", N_I * N_TS);

    let state = Arc::new(AppState { flux_luts, thermal_luts });

    let app = Router::new()
        .route("/echo", get(echo_handler))
        .with_state(state);

    let addr = "0.0.0.0:1416";
    let listener = TcpListener::bind(addr).await?;
    tracing::info!("Optec Thermal API listening on {addr}");
    axum::serve(listener, app).await?;

    Ok(())
}
