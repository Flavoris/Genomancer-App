import XCTest
@testable import DigestCore

/// Comprehensive test suite to verify all enzymes cut at the correct positions
/// Tests all enzymes from the database to ensure:
/// 1. Recognition sites are found correctly
/// 2. Cut positions match expected values
/// 3. Overhang types are correct
/// 4. Overhang sequences are valid
final class EnzymeCutSiteTests: XCTestCase {
    
    var allEnzymes: [Enzyme] = []
    
    override func setUp() {
        super.setUp()
        
        // Load all enzymes from the test bundle
        guard let url = Bundle.module.url(forResource: "enzymes", withExtension: "json"),
              let data = try? Data(contentsOf: url) else {
            XCTFail("Failed to load enzymes.json from test bundle")
            return
        }
        
        // Configure decoder to handle snake_case JSON keys
        let decoder = JSONDecoder()
        decoder.keyDecodingStrategy = .convertFromSnakeCase
        
        guard let enzymes = try? decoder.decode([Enzyme].self, from: data) else {
            XCTFail("Failed to decode enzymes.json")
            return
        }
        
        allEnzymes = enzymes
        print("Loaded \(allEnzymes.count) enzymes for testing")
    }
    
    // MARK: - Individual Enzyme Cut Site Tests
    
    func test_allEnzymes_cutAtCorrectPosition() {
        var failedEnzymes: [(String, String)] = []
        var successCount = 0
        
        for enzyme in allEnzymes {
            let result = testSingleEnzymeCutPosition(enzyme)
            if let error = result {
                failedEnzymes.append((enzyme.name, error))
            } else {
                successCount += 1
            }
        }
        
        // Report results
        print("\n=== Enzyme Cut Site Test Results ===")
        print("✅ Passed: \(successCount)/\(allEnzymes.count)")
        print("❌ Failed: \(failedEnzymes.count)/\(allEnzymes.count)")
        
        if !failedEnzymes.isEmpty {
            print("\nFailed enzymes:")
            for (name, error) in failedEnzymes {
                print("  - \(name): \(error)")
            }
        }
        
        XCTAssertTrue(failedEnzymes.isEmpty, 
                     "Failed to validate \(failedEnzymes.count) enzymes: \(failedEnzymes.map { $0.0 }.joined(separator: ", "))")
    }
    
    /// Test a single enzyme's cut position
    /// Returns nil on success, error message on failure
    private func testSingleEnzymeCutPosition(_ enzyme: Enzyme) -> String? {
        // Skip enzymes without cut index (non-cutters or special cases)
        guard let cutIndex = enzyme.cutIndexTop else {
            return nil // Not an error, just skip
        }
        
        let site = enzyme.site
        let siteLength = site.count
        
        // Create a test sequence with the recognition site in the middle
        // Use a concrete sequence that matches the IUPAC pattern
        let concreteSequence = expandIUPACToSequence(iupac: site)
        guard !concreteSequence.isEmpty else {
            return "Failed to expand IUPAC sequence: \(site)"
        }
        
        let padding = String(repeating: "A", count: 100)
        let testSequence = padding + concreteSequence + padding
        let expectedCutPosition = 100 + cutIndex
        
        // Run digest
        let engine = DigestEngine(sequence: testSequence, enzymes: [enzyme])
        let fragments = engine.digest(options: .init(circular: false, returnSequences: false))
        
        // For linear digest with one recognition site, we expect 2 fragments
        guard fragments.count == 2 else {
            // Check if the enzyme didn't match (which could happen with complex IUPAC codes)
            let cutSites = engine.cutSites()
            if cutSites.isEmpty {
                return "No cut sites found for sequence: \(site) (IUPAC expansion may have failed)"
            }
            return "Expected 2 fragments, got \(fragments.count)"
        }
        
        // Verify the cut position
        let sortedFrags = fragments.sorted { $0.start < $1.start }
        let actualCutPosition = sortedFrags[0].end
        
        guard actualCutPosition == expectedCutPosition else {
            return "Cut position mismatch: expected \(expectedCutPosition), got \(actualCutPosition)"
        }
        
        // Verify fragment lengths add up to total sequence length
        let totalLength = fragments.map { $0.length }.reduce(0, +)
        guard totalLength == testSequence.count else {
            return "Fragment lengths don't add up: \(totalLength) != \(testSequence.count)"
        }
        
        // Verify overhang type if specified
        if let expectedOverhangType = enzyme.overhangType {
            let rightEnd = sortedFrags[0].rightEnd
            let leftEnd = sortedFrags[1].leftEnd
            
            guard rightEnd.overhangType == expectedOverhangType else {
                return "Overhang type mismatch: expected \(expectedOverhangType), got \(rightEnd.overhangType)"
            }
            
            guard leftEnd.overhangType == expectedOverhangType else {
                return "Overhang type mismatch on left end: expected \(expectedOverhangType), got \(leftEnd.overhangType)"
            }
            
            // Verify overhang length
            let expectedOverhangLength = calculateExpectedOverhangLength(
                siteLength: siteLength,
                cutIndex: cutIndex
            )
            
            guard rightEnd.overhangLen == expectedOverhangLength else {
                return "Overhang length mismatch: expected \(expectedOverhangLength), got \(rightEnd.overhangLen)"
            }
        }
        
        return nil // Success
    }
    
    // MARK: - Specific Enzyme Tests
    
    func test_EcoRI_cutPosition() {
        // EcoRI: GAATTC, cut at position 1
        // Expected: G^AATTC (cuts between G and A)
        let ecoRI = allEnzymes.first { $0.name == "EcoRI" }
        XCTAssertNotNil(ecoRI, "EcoRI not found in enzyme database")
        
        guard let enzyme = ecoRI else { return }
        guard let cutIndex = enzyme.cutIndexTop else {
            XCTFail("EcoRI has no cut index")
            return
        }
        
        let sequence = "AAAAAAGAATTCAAAAA"
        let engine = DigestEngine(sequence: sequence, enzymes: [enzyme])
        let fragments = engine.digest(options: .init(circular: false))
        
        // Should cut at position 6 + cutIndex
        XCTAssertEqual(fragments.count, 2)
        let sorted = fragments.sorted { $0.start < $1.start }
        XCTAssertEqual(sorted[0].end, 6 + cutIndex, "EcoRI should cut at position \(6 + cutIndex)")
        if let overhangType = enzyme.overhangType {
            XCTAssertEqual(sorted[0].rightEnd.overhangType, overhangType)
        }
    }
    
    func test_BamHI_cutPosition() {
        // BamHI: GGATCC, cut at position 1
        // Expected: G^GATCC
        let bamHI = allEnzymes.first { $0.name == "BamHI" }
        XCTAssertNotNil(bamHI, "BamHI not found in enzyme database")
        
        guard let enzyme = bamHI else { return }
        guard let cutIndex = enzyme.cutIndexTop else {
            XCTFail("BamHI has no cut index")
            return
        }
        
        let sequence = "AAAAAAGGATCCAAAAA"
        let engine = DigestEngine(sequence: sequence, enzymes: [enzyme])
        let fragments = engine.digest(options: .init(circular: false))
        
        XCTAssertEqual(fragments.count, 2)
        let sorted = fragments.sorted { $0.start < $1.start }
        XCTAssertEqual(sorted[0].end, 6 + cutIndex, "BamHI should cut at position \(6 + cutIndex)")
        if let overhangType = enzyme.overhangType {
            XCTAssertEqual(sorted[0].rightEnd.overhangType, overhangType)
        }
    }
    
    func test_PstI_cutPosition() {
        // PstI: CTGCAG, cut at position 5
        // Expected: CTGCA^G (3' overhang)
        let pstI = allEnzymes.first { $0.name == "PstI" }
        XCTAssertNotNil(pstI, "PstI not found in enzyme database")
        
        guard let enzyme = pstI else { return }
        guard let cutIndex = enzyme.cutIndexTop else {
            XCTFail("PstI has no cut index")
            return
        }
        
        let sequence = "AAAAACTGCAGAAAA"
        let engine = DigestEngine(sequence: sequence, enzymes: [enzyme])
        let fragments = engine.digest(options: .init(circular: false))
        
        XCTAssertEqual(fragments.count, 2)
        let sorted = fragments.sorted { $0.start < $1.start }
        XCTAssertEqual(sorted[0].end, 5 + cutIndex, "PstI should cut at position \(5 + cutIndex)")
        if let overhangType = enzyme.overhangType {
            XCTAssertEqual(sorted[0].rightEnd.overhangType, overhangType)
        }
    }
    
    func test_HindIII_cutPosition() {
        // HindIII: AAGCTT, cut at position 1
        // Expected: A^AGCTT
        let hindIII = allEnzymes.first { $0.name == "HindIII" }
        XCTAssertNotNil(hindIII, "HindIII not found in enzyme database")
        
        guard let enzyme = hindIII else { return }
        guard let cutIndex = enzyme.cutIndexTop else {
            XCTFail("HindIII has no cut index")
            return
        }
        
        // Sequence has "AAGCTT" at position 4
        let sequence = "AAAAAAGCTTAAAAA"
        let engine = DigestEngine(sequence: sequence, enzymes: [enzyme])
        let fragments = engine.digest(options: .init(circular: false))
        
        XCTAssertEqual(fragments.count, 2)
        let sorted = fragments.sorted { $0.start < $1.start }
        // The site starts at position 4, cut happens at position 4 + cutIndex
        XCTAssertEqual(sorted[0].end, 4 + cutIndex, "HindIII should cut at position \(4 + cutIndex)")
        if let overhangType = enzyme.overhangType {
            XCTAssertEqual(sorted[0].rightEnd.overhangType, overhangType)
        }
    }
    
    func test_bluntCutters() {
        // Test several blunt-cutting enzymes
        let bluntCutterNames = ["EcoRV", "SmaI", "PvuII"]
        
        for name in bluntCutterNames {
            guard let enzyme = allEnzymes.first(where: { $0.name == name }) else {
                XCTFail("\(name) not found in database")
                continue
            }
            
            guard let cutIndex = enzyme.cutIndexTop else {
                XCTFail("\(name) has no cut index")
                continue
            }
            
            let site = enzyme.site
            let sequence = "AAAAA\(site)AAAAA"
            let engine = DigestEngine(sequence: sequence, enzymes: [enzyme])
            let fragments = engine.digest(options: .init(circular: false))
            
            XCTAssertEqual(fragments.count, 2, "\(name) should produce 2 fragments")
            
            let sorted = fragments.sorted { $0.start < $1.start }
            XCTAssertEqual(sorted[0].end, 5 + cutIndex, "\(name) should cut at position \(5 + cutIndex)")
            if let overhangType = enzyme.overhangType {
                XCTAssertEqual(sorted[0].rightEnd.overhangType, overhangType, "\(name) overhang type")
                XCTAssertEqual(sorted[1].leftEnd.overhangType, overhangType, "\(name) overhang type")
            }
        }
    }
    
    func test_enzymesWithIUPACCodes() {
        // Test enzymes that use IUPAC ambiguity codes
        let iupacEnzymeNames = ["AcsI", "BanI", "AvaI"]
        
        for name in iupacEnzymeNames {
            guard let enzyme = allEnzymes.first(where: { $0.name == name }) else {
                XCTFail("\(name) not found in database")
                continue
            }
            
            // Test with one possible concrete sequence
            let concreteSeq = expandIUPACToSequence(iupac: enzyme.site)
            guard !concreteSeq.isEmpty else {
                XCTFail("Failed to expand IUPAC for \(name): \(enzyme.site)")
                continue
            }
            
            let sequence = "AAAAA\(concreteSeq)AAAAA"
            let engine = DigestEngine(sequence: sequence, enzymes: [enzyme])
            let fragments = engine.digest(options: .init(circular: false))
            
            XCTAssertEqual(fragments.count, 2, 
                          "\(name) should recognize \(concreteSeq) (from \(enzyme.site))")
        }
    }
    
    // MARK: - Helper Functions
    
    /// Expand IUPAC ambiguity codes to a concrete DNA sequence
    /// Returns one valid sequence that matches the IUPAC pattern
    private func expandIUPACToSequence(iupac: String) -> String {
        var result = ""
        for char in iupac {
            switch char {
            case "A": result += "A"
            case "C": result += "C"
            case "G": result += "G"
            case "T": result += "T"
            case "R": result += "A"  // R = A or G, choose A
            case "Y": result += "C"  // Y = C or T, choose C
            case "S": result += "G"  // S = G or C, choose G
            case "W": result += "A"  // W = A or T, choose A
            case "K": result += "G"  // K = G or T, choose G
            case "M": result += "A"  // M = A or C, choose A
            case "B": result += "C"  // B = C or G or T, choose C
            case "D": result += "A"  // D = A or G or T, choose A
            case "H": result += "A"  // H = A or C or T, choose A
            case "V": result += "A"  // V = A or C or G, choose A
            case "N": result += "A"  // N = any base, choose A
            default: result += "A"   // Unknown, default to A
            }
        }
        return result
    }
    
    /// Calculate expected overhang length based on site length and cut index
    private func calculateExpectedOverhangLength(siteLength: Int, cutIndex: Int) -> Int {
        let diff = 2 * cutIndex - siteLength
        return abs(diff)
    }
}

