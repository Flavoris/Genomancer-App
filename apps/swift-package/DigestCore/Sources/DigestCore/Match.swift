import Foundation

/// Optimized version: works with pre-uppercased character arrays (no allocations)
public func findAllMatches(sequenceChars: [Character], motifChars: [Character]) -> [Int] {
    guard sequenceChars.count >= motifChars.count, motifChars.count > 0 else { return [] }

    var hits: [Int] = []
    outer: for i in 0...(sequenceChars.count - motifChars.count) {
        for j in 0..<motifChars.count {
            if !IUPAC.matchesFast(iupac: motifChars[j], base: sequenceChars[i+j]) { continue outer }
        }
        hits.append(i)
    }
    return hits
}

/// Legacy version: for backward compatibility or when strings are not pre-processed
public func findAllMatches(sequence: String, motifIUPAC: String) -> [Int] {
    let seq = Array(sequence.uppercased())
    let motif = Array(motifIUPAC.uppercased())
    return findAllMatches(sequenceChars: seq, motifChars: motif)
}

public func findMatchesCircular(sequence: String, motifIUPAC: String) -> [Int] {
    let n = sequence.count
    guard n > 0 else { return [] }
    let ext = sequence + sequence.prefix(max(0, motifIUPAC.count - 1))
    return findAllMatches(sequence: String(ext), motifIUPAC: motifIUPAC)
        .filter { $0 < n }
}

