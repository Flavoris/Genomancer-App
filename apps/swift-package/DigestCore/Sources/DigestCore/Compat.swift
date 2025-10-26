import Foundation

// ============================================================================
// MARK: - Data Structures
// ============================================================================

/// Information about a fragment end with ligation compatibility data
public struct FragmentEndInfo: Hashable {
    public let fragmentId: Int
    public let polarity: String  // "left" or "right"
    public let position: Int
    public let endInfo: EndInfo
    
    public init(fragmentId: Int, polarity: String, position: Int, endInfo: EndInfo) {
        self.fragmentId = fragmentId
        self.polarity = polarity
        self.position = position
        self.endInfo = endInfo
    }
}

/// Result of compatibility check between two ends
public struct CompatibilityResult: Hashable, Identifiable {
    public let id = UUID()
    public let endA: FragmentEndInfo
    public let endB: FragmentEndInfo
    public let compatible: Bool
    public let directional: Bool
    public let gcPercentA: Double
    public let gcPercentB: Double
    public let tmA: Double
    public let tmB: Double
    public let note: String
    
    public init(endA: FragmentEndInfo, endB: FragmentEndInfo, compatible: Bool,
                directional: Bool, gcPercentA: Double, gcPercentB: Double,
                tmA: Double, tmB: Double, note: String) {
        self.endA = endA
        self.endB = endB
        self.compatible = compatible
        self.directional = directional
        self.gcPercentA = gcPercentA
        self.gcPercentB = gcPercentB
        self.tmA = tmA
        self.tmB = tmB
        self.note = note
    }
    
    public func hash(into hasher: inout Hasher) {
        hasher.combine(endA)
        hasher.combine(endB)
        hasher.combine(compatible)
        hasher.combine(directional)
        hasher.combine(gcPercentA)
        hasher.combine(gcPercentB)
        hasher.combine(tmA)
        hasher.combine(tmB)
        hasher.combine(note)
    }
    
    public static func == (lhs: CompatibilityResult, rhs: CompatibilityResult) -> Bool {
        lhs.endA == rhs.endA &&
        lhs.endB == rhs.endB &&
        lhs.compatible == rhs.compatible &&
        lhs.directional == rhs.directional &&
        lhs.gcPercentA == rhs.gcPercentA &&
        lhs.gcPercentB == rhs.gcPercentB &&
        lhs.tmA == rhs.tmA &&
        lhs.tmB == rhs.tmB &&
        lhs.note == rhs.note
    }
    
    /// Calculate compatibility strength score (0.0 to 1.0)
    /// Based on GC content, Tm, overhang length, and directionality
    public var strengthScore: Double {
        var score = 0.0
        let weights = (gc: 0.2, tm: 0.3, length: 0.3, directional: 0.2)
        
        // GC content score (optimal around 50%)
        let avgGC = (gcPercentA + gcPercentB) / 2.0
        let gcScore = 1.0 - abs(avgGC - 50.0) / 50.0
        score += gcScore * weights.gc
        
        // Tm score (higher is better, normalize to 0-1 with 50°C as max)
        let avgTm = (tmA + tmB) / 2.0
        let tmScore = min(avgTm / 50.0, 1.0)
        score += tmScore * weights.tm
        
        // Length score (longer is better, normalize with 6bp as ideal)
        let length = Double(endA.endInfo.overhangLen)
        let lengthScore = min(length / 6.0, 1.0)
        score += lengthScore * weights.length
        
        // Directionality bonus
        if directional {
            score += weights.directional
        }
        
        return min(max(score, 0.0), 1.0)
    }
    
    /// Color representing the strength of this compatibility
    public var strengthColor: (red: Double, green: Double, blue: Double, alpha: Double) {
        if strengthScore >= 0.75 {
            return (0.0, 0.8, 0.0, 1.0)  // Green
        } else if strengthScore >= 0.5 {
            return (0.0, 0.8, 1.0, 1.0)  // Genomancer accent (cyan-ish)
        } else {
            return (1.0, 0.6, 0.0, 1.0)  // Orange
        }
    }
    
    /// Label describing the strength of this compatibility
    public var strengthLabel: String {
        if strengthScore >= 0.75 {
            return "Excellent"
        } else if strengthScore >= 0.5 {
            return "Good"
        } else if strengthScore >= 0.25 {
            return "Moderate"
        } else {
            return "Weak"
        }
    }
}

// ============================================================================
// MARK: - Core Compatibility Functions
// ============================================================================

/// Determines if two fragment ends are compatible for ligation based on concrete sequence matching
/// - Parameters:
///   - a: First fragment end
///   - b: Second fragment end
///   - includeBlunt: Whether to consider blunt-blunt ends as compatible
///   - minOverhang: Minimum overhang length required for sticky ends
/// - Returns: True if the ends can ligate together
public func compatibleConcrete(a: EndInfo, b: EndInfo,
                               includeBlunt: Bool = false,
                               minOverhang: Int = 1) -> Bool {
    // sticky↔blunt always false
    let aSticky = (a.overhangLen > 0 && a.overhangSeq5to3?.isEmpty == false)
    let bSticky = (b.overhangLen > 0 && b.overhangSeq5to3?.isEmpty == false)
    if aSticky != bSticky { return false }
    if !aSticky && !bSticky { return includeBlunt } // blunt↔blunt
    
    if a.overhangType != b.overhangType { return false }
    if a.overhangLen < minOverhang || b.overhangLen < minOverhang { return false }
    if a.overhangLen != b.overhangLen { return false }
    
    guard let sa = a.overhangSeq5to3, let sb = b.overhangSeq5to3 else { return false }
    return sa.uppercased() == revcomp(sb).uppercased()
}

/// Check if a compatible pair enforces directionality
/// - Parameters:
///   - a: First fragment end
///   - b: Second fragment end
/// - Returns: True if directional (non-palindromic)
public func isDirectional(a: EndInfo, b: EndInfo) -> Bool {
    guard let s = a.overhangSeq5to3 else { return false }
    let upper = s.uppercased()
    return upper != revcomp(upper)
}

// ============================================================================
// MARK: - Heuristic Calculations
// ============================================================================

/// Calculate GC percentage of a DNA sequence
/// - Parameter seq: DNA sequence
/// - Returns: GC percentage (0-100)
public func calculateGCPercent(_ seq: String) -> Double {
    guard !seq.isEmpty else { return 0.0 }
    let upper = seq.uppercased()
    let gcCount = upper.filter { $0 == "G" || $0 == "C" }.count
    return (Double(gcCount) / Double(seq.count)) * 100.0
}

/// Calculate rough melting temperature using Wallace rule
/// Tm ≈ 2*(A+T) + 4*(G+C)
/// This is only accurate for short oligonucleotides (<14 nt)
/// - Parameter seq: DNA sequence
/// - Returns: Estimated Tm in °C
public func calculateTm(_ seq: String) -> Double {
    guard !seq.isEmpty else { return 0.0 }
    let upper = seq.uppercased()
    let atCount = upper.filter { $0 == "A" || $0 == "T" }.count
    let gcCount = upper.filter { $0 == "G" || $0 == "C" }.count
    return 2.0 * Double(atCount) + 4.0 * Double(gcCount)
}

// ============================================================================
// MARK: - Comprehensive Compatibility Analysis
// ============================================================================

/// Calculate compatibility between all pairs of fragment ends
/// - Parameters:
///   - fragments: Array of fragments to analyze
///   - includeBlunt: Include blunt-blunt compatibility
///   - minOverhang: Minimum overhang length
///   - requireDirectional: Only return directional pairs
/// - Returns: Array of CompatibilityResult objects
public func calculateCompatibility(fragments: [Fragment],
                                   includeBlunt: Bool = false,
                                   minOverhang: Int = 1,
                                   requireDirectional: Bool = false) -> [CompatibilityResult] {
    // Build list of all fragment ends
    var ends: [FragmentEndInfo] = []
    for (index, fragment) in fragments.enumerated() {
        // Left end
        ends.append(FragmentEndInfo(
            fragmentId: index,
            polarity: "left",
            position: fragment.start,
            endInfo: fragment.leftEnd
        ))
        // Right end
        ends.append(FragmentEndInfo(
            fragmentId: index,
            polarity: "right",
            position: fragment.end,
            endInfo: fragment.rightEnd
        ))
    }
    
    var results: [CompatibilityResult] = []
    
    // Check all pairs
    for i in 0..<ends.count {
        for j in (i+1)..<ends.count {
            let endA = ends[i]
            let endB = ends[j]
            
            // Check compatibility
            let compatible = compatibleConcrete(
                a: endA.endInfo,
                b: endB.endInfo,
                includeBlunt: includeBlunt,
                minOverhang: minOverhang
            )
            
            guard compatible else { continue }
            
            // Check directionality
            let directional = isDirectional(a: endA.endInfo, b: endB.endInfo)
            
            // Skip if directionality is required but not met
            if requireDirectional && !directional {
                continue
            }
            
            // Calculate heuristics
            let seqA = endA.endInfo.overhangSeq5to3 ?? ""
            let seqB = endB.endInfo.overhangSeq5to3 ?? ""
            let gcA = calculateGCPercent(seqA)
            let gcB = calculateGCPercent(seqB)
            let tmA = calculateTm(seqA)
            let tmB = calculateTm(seqB)
            
            // Generate note
            let note = generateCompatibilityNote(endA: endA.endInfo, endB: endB.endInfo, directional: directional)
            
            results.append(CompatibilityResult(
                endA: endA,
                endB: endB,
                compatible: true,
                directional: directional,
                gcPercentA: gcA,
                gcPercentB: gcB,
                tmA: tmA,
                tmB: tmB,
                note: note
            ))
        }
    }
    
    return results
}

/// Generate a descriptive note about the compatibility
/// - Parameters:
///   - endA: First end info
///   - endB: Second end info
///   - directional: Whether the pair is directional
/// - Returns: Descriptive note
public func generateCompatibilityNote(endA: EndInfo, endB: EndInfo, directional: Bool) -> String {
    if endA.overhangLen == 0 {
        return "Blunt-blunt ligation"
    }
    
    let dirNote = directional ? "directional" : "non-directional (palindromic)"
    let typeStr = endA.overhangType == .fivePrime ? "5' overhang" : "3' overhang"
    return "\(typeStr), \(endA.overhangLen) bp overhang, \(dirNote)"
}

