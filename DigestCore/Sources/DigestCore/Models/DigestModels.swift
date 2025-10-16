import Foundation

public enum Strand {
    case top
    case bottom
}

public struct CutSite: Hashable {
    public let enzyme: Enzyme
    public let position: Int              // absolute 0-based
    public let strand: Strand
    
    public init(enzyme: Enzyme, position: Int, strand: Strand) {
        self.enzyme = enzyme
        self.position = position
        self.strand = strand
    }
}

public struct EndInfo: Hashable {
    public let overhangType: OverhangType
    public let overhangSeq5to3: String?
    public let sourceEnzyme: String?
    
    public init(overhangType: OverhangType, overhangSeq5to3: String?, sourceEnzyme: String?) {
        self.overhangType = overhangType
        self.overhangSeq5to3 = overhangSeq5to3
        self.sourceEnzyme = sourceEnzyme
    }
}

public struct Fragment: Hashable {
    public let start: Int
    public let end: Int                   // half-open [start, end)
    public let length: Int
    public let leftEnd: EndInfo
    public let rightEnd: EndInfo
    public let sequence: String?
    
    public init(
        start: Int,
        end: Int,
        length: Int,
        leftEnd: EndInfo,
        rightEnd: EndInfo,
        sequence: String?
    ) {
        self.start = start
        self.end = end
        self.length = length
        self.leftEnd = leftEnd
        self.rightEnd = rightEnd
        self.sequence = sequence
    }
}

