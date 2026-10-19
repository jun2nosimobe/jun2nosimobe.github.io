// ecm_worker.js

// Wasmが読み込まれた瞬間にメインスレッドに通知する設定
self.Module = {
    onRuntimeInitialized: function() {
        postMessage({ type: 'ready' });
    }
};

// C++をコンパイルしたWasmモジュールを読み込む
importScripts('ecm_module.js');

onmessage = function(e) {
    if (e.data.cmd === 'run_ecm') {
        const { n_str, B1, seed, num_curves, workerId } = e.data;
        
        // C++のバッチ処理を実行
        const result = Module.run_ecm_batch(n_str, B1, seed, num_curves);
        
        // メインスレッドに結果を返却
        postMessage({
            type: 'result',
            workerId: workerId,
            status: result.status,
            factor: result.factor,
            phase: result.phase,
            prime: result.prime,
            k: result.k,
            Z: result.Z,
            A: result.A,
            B: result.B,
            P_x: result.P_x,
            P_y: result.P_y,
            curves_tried: result.curves_tried
        });
    } else if (e.data.cmd === 'run_rho') {
        // ロー法 (メインスレッドを固めないようWorkerで実行)
        const { n_str, steps, workerId } = e.data;
        postMessage({
            type: 'rho_result',
            workerId: workerId,
            factor: Module.pollard_rho(n_str, steps)
        });
    }
};