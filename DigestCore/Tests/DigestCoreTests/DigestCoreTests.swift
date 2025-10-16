import XCTest
@testable import DigestCore

final class DigestCoreTests: XCTestCase {
    
    func testDigestCoreInitialization() throws {
        let digestCore = DigestCore()
        XCTAssertNotNil(digestCore)
    }
    
    func testAnalyzeSequence() throws {
        let digestCore = DigestCore()
        let result = digestCore.analyzeSequence("ATCG")
        XCTAssertEqual(result, "Analyzing sequence: ATCG")
    }
    
    func testFindAllMatches() throws {
        // Test exact match
        let matches = findAllMatches(sequence: "ATGAATTCGC", motifIUPAC: "GAATTC")
        XCTAssertEqual(matches, [2], "Should find EcoRI site at position 2")
        
        // Test IUPAC degenerate match (R = A or G)
        let degenerateMatches = findAllMatches(sequence: "ATGAATTC", motifIUPAC: "RTG")
        XCTAssertEqual(degenerateMatches, [0], "Should match ATG with RTG")
    }
    
    func testFindMatchesCircular() throws {
        // Test circular match wrapping around
        let circularMatches = findMatchesCircular(sequence: "TCATGA", motifIUPAC: "GATC")
        XCTAssertEqual(circularMatches, [4], "Should find GATC wrapping from position 4")
    }
}
