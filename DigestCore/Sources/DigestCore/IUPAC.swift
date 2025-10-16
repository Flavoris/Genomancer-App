import Foundation

public enum IUPAC {
    public static let map: [Character: Set<Character>] = [
        "A":["A"], "C":["C"], "G":["G"], "T":["T"],
        "R":["A","G"], "Y":["C","T"], "W":["A","T"], "S":["G","C"],
        "M":["A","C"], "K":["G","T"], "B":["C","G","T"], "D":["A","G","T"],
        "H":["A","C","T"], "V":["A","C","G"], "N":["A","C","G","T"]
    ]
    public static func matches(iupac: Character, base: Character) -> Bool {
        guard let set = map[Character(iupac.uppercased())] else { return false }
        return set.contains(Character(base.uppercased()))
    }
}

public func rc(_ s: String) -> String {
    let comp: [Character:Character] = ["A":"T","T":"A","C":"G","G":"C",
                                       "a":"t","t":"a","c":"g","g":"c"]
    return String(s.reversed().map{ comp[$0] ?? $0 })
}

/// Reverse complement with full IUPAC support
public func revcomp(_ s: String) -> String {
    let comp: [Character:Character] = ["A":"T","T":"A","C":"G","G":"C",
                                       "a":"t","t":"a","c":"g","g":"c",
                                       "R":"Y","Y":"R","S":"S","W":"W",
                                       "K":"M","M":"K","B":"V","D":"H",
                                       "H":"D","V":"B","N":"N",
                                       "r":"y","y":"r","s":"s","w":"w","k":"m","m":"k",
                                       "b":"v","d":"h","h":"d","v":"b","n":"n"]
    return String(s.reversed().map { comp[$0] ?? $0 })
}

/// IUPAC compatibility: at each column, A's set ∩ complement(B's set) ≠ ∅
public func iupacCompatible(a: String, b: String) -> Bool {
    let A = Array(a.uppercased()), B = Array(b.uppercased())
    guard A.count == B.count else { return false }
    for i in 0..<A.count {
        guard let SA = IUPAC.map[A[i]], let SB = IUPAC.map[B[i]] else { return false }
        let SBc = Set(SB.map { revcomp(String($0)).first! }) // G→C, etc.
        if SA.intersection(SBc).isEmpty { return false }
    }
    return true
}

