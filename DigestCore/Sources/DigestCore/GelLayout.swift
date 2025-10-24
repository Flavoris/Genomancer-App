// GelLayout.swift
import Foundation
import CoreGraphics

/// Physics-ish gel mobility utilities using a Ferguson model.
/// log(mu(C)) = log(mu0(bp)) - K(bp) * C
/// We normalize mobility across [minBP, maxBP] to an Rf in [0,1] and map to pixels.
public enum GelLayout {

    // ---- Tunables (visual → calibrate later with real ladders) ----
    // Zero-gel baseline mobility ~ A * bp^(-alpha)
    private static let A: Double = 1.10
    private static let alpha: Double = 0.40   // slightly flatter size-dependence

    // Retardation coefficient K(bp) ~ k0 + k1 * log10(bp)
    private static let k0: Double = -0.02
    private static let k1: Double =  0.18    // gentler % slow-down growth with size

    // MARK: - Ferguson helpers

    public static func mu0(bp: Int) -> Double {
        let b = max(1, bp)
        return A * pow(Double(b), -alpha)
    }

    public static func K(bp: Int) -> Double {
        let b = max(1, bp)
        return k0 + k1 * log10(Double(b))
    }

    private static func sizeFactor(bp: Int) -> Double {
        let b = max(2, bp) // avoid log(1)
        // Smaller fragments should move faster → larger factor.
        // 1 / log10(bp) grows as bp shrinks; clamp to avoid extremes.
        let f = 1.0 / log10(Double(b))
        return min(2.0, max(0.5, f))
    }

    /// Mobility at gel percentage C (e.g., 1.0, 2.0, 3.0).
    /// Fragment size contributes both via μ₀/K (Ferguson model) and an explicit size factor
    /// for clearer, tunable size-dependence. Smaller fragments move faster.
    public static func mobility(bp: Int, gelPercent C: Double) -> Double {
        let ferguson = mu0(bp: bp) * exp(-K(bp: bp) * max(0.0, C))
        return ferguson * sizeFactor(bp: bp)
    }

    /// Map bp → y using Ferguson mobility in log-space, normalized against a reference gel%.
    /// Bands actually migrate when gel% changes, anchored to a 2.0% reference window.
    /// Log-space spreads out near-zero mobilities for large fragments, preventing quantization at the top.
    /// Larger fragments (bp↑) end up nearer the top (y closer to `top`).
    public static func y(bp: Int, minBP: Int, maxBP: Int,
                         top: CGFloat, height: CGFloat, gelPercent: Double) -> CGFloat
    {
        let lo = min(minBP, maxBP)
        let hi = max(minBP, maxBP)

        let step = max(1, (hi - lo) / 18)
        let eps  = 1e-12

        // --- Reference window (fixed gel%) for normalization ---
        let Cref = 2.0  // pick the midpoint of your UI presets
        let logsRef = stride(from: lo, through: hi, by: step).map {
            log( max(eps, mobility(bp: $0, gelPercent: Cref)) )
        }
        guard let lMinRef = logsRef.min(), let lMaxRef = logsRef.max(), lMaxRef > lMinRef else {
            return top + height
        }
        let denom = lMaxRef - lMinRef

        // Current band position uses the CURRENT gel%
        let lmu = log( max(eps, mobility(bp: bp, gelPercent: gelPercent)) )
        var Rf  = (lmu - lMinRef) / denom
        // Clamp: if outside due to extreme C, keep visible
        Rf = min(1, max(0, Rf))
        return top + CGFloat(Rf) * height
    }
}
