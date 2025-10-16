import XCTest
@testable import DigestCore

final class DigestCoreTests: XCTestCase {
    func test_IUPAC_match() {
        XCTAssertTrue(IUPAC.matches(iupac: "R", base: "A"))
        XCTAssertFalse(IUPAC.matches(iupac: "R", base: "C"))
    }

    func test_EcoRI_linear_digest() {
        let ecoRI = Enzyme(name: "EcoRI", site: "GAATTC", cutIndexTop: 1, cutIndexBottom: 5, overhangType: .fivePrime, notes: nil)
        let engine = DigestEngine(sequence: "AAGAATTCTT", enzymes: [ecoRI])
        let frags = engine.digest(options: .init(circular: false))
        XCTAssertEqual(frags.map(\.length).sorted(), [3, 7])
    }

    func test_EcoRI_overhangs_and_compat() {
        // EcoRI: GAATTC with top cut at 1, bottom 5 → 5' overhang length 4 (AATT)
        let ecoRI = Enzyme(name: "EcoRI", site: "GAATTC",
                           cutIndexTop: 1, cutIndexBottom: 5,
                           overhangType: .fivePrime, notes: nil)

        let seq = "AAGAATTCTT" // ...GAATTC...
        let engine = DigestEngine(sequence: seq, enzymes: [ecoRI])
        let frags = engine.digest(options: .init(circular: false, returnSequences: false))
        // Expect 2 fragments from single cut at position 3
        XCTAssertEqual(frags.map(\.length).sorted(), [3, 7])

        // Find the fragment whose left or right end was created by EcoRI cut
        let ends = frags.flatMap { [$0.leftEnd, $0.rightEnd] }
        // Sticky ends should have AATT, len 4, 5' type
        let stickies = ends.filter { $0.overhangLen > 0 }
        XCTAssertFalse(stickies.isEmpty)
        for e in stickies {
            XCTAssertEqual(e.overhangType, .fivePrime)
            XCTAssertEqual(e.overhangLen, 4)
            XCTAssertEqual(e.overhangSeq5to3, "AATT")
        }
        // Compatibility: two EcoRI ends should be mutually compatible (revcomp(AATT) == AATT' on the other end)
        XCTAssertTrue(stickies.count >= 2)
        XCTAssertTrue(compatibleConcrete(a: stickies[0], b: stickies[1]))
    }

    func test_PstI_overhangs() {
        // PstI: CTGCAG with top cut at 5, bottom at 1 → 3' overhang length 4 ("TGCA")
        let pstI = Enzyme(name: "PstI", site: "CTGCAG",
                          cutIndexTop: 5, cutIndexBottom: 1,
                          overhangType: .threePrime, notes: nil)
        let seq = "AAACTGCAGTT" // ...CTGCAG...
        let engine = DigestEngine(sequence: seq, enzymes: [pstI])
        let frags = engine.digest(options: .init(circular: false))
        let ends = frags.flatMap { [$0.leftEnd, $0.rightEnd] }
        let stickies = ends.filter { $0.overhangLen > 0 }
        XCTAssertFalse(stickies.isEmpty)
        for e in stickies {
            XCTAssertEqual(e.overhangType, .threePrime)
            XCTAssertEqual(e.overhangLen, 4)
            XCTAssertEqual(e.overhangSeq5to3, "TGCA")
        }
    }

    // MARK: - Task 9.7: QA Regression Tests
    
    func test_mixedDigest_EcoRI_PstI_2kb() {
        // Create a ~2 kb toy sequence with known EcoRI and PstI sites
        // EcoRI: GAATTC (cut at 1 → 5' overhang AATT)
        // PstI: CTGCAG (cut at 5 → 3' overhang TGCA)
        
        let ecoRI = Enzyme(name: "EcoRI", site: "GAATTC",
                           cutIndexTop: 1, cutIndexBottom: 5,
                           overhangType: .fivePrime, notes: nil)
        let pstI = Enzyme(name: "PstI", site: "CTGCAG",
                          cutIndexTop: 5, cutIndexBottom: 1,
                          overhangType: .threePrime, notes: nil)
        
        // Build a 2kb sequence with specific sites at known positions
        // Position 500: EcoRI site (GAATTC)
        // Position 1200: PstI site (CTGCAG)
        // Position 1800: EcoRI site (GAATTC)
        var seq = String(repeating: "A", count: 500)
        seq += "GAATTC"  // EcoRI at 500, length 6
        seq += String(repeating: "T", count: 688)  // 500 + 6 + 688 = 1194
        seq += "CTGCAG"  // PstI at 1194, length 6
        seq += String(repeating: "G", count: 594)  // 1194 + 6 + 594 = 1794
        seq += "GAATTC"  // EcoRI at 1794, length 6
        seq += String(repeating: "C", count: 200)  // 1794 + 6 + 200 = 2000
        
        XCTAssertEqual(seq.count, 2000, "Test sequence should be 2kb")
        
        let engine = DigestEngine(sequence: seq, enzymes: [ecoRI, pstI])
        let frags = engine.digest(options: .init(circular: false, returnSequences: true))
        
        // Expected cuts:
        // - EcoRI at position 500 + 1 = 501
        // - PstI at position 1194 + 5 = 1199
        // - EcoRI at position 1794 + 1 = 1795
        // Expected fragments:
        // 1. [0, 501) = 501 bp
        // 2. [501, 1199) = 698 bp
        // 3. [1199, 1795) = 596 bp
        // 4. [1795, 2000) = 205 bp
        
        XCTAssertEqual(frags.count, 4, "Should produce 4 fragments")
        
        let sortedFrags = frags.sorted { $0.start < $1.start }
        
        // Fragment 1: natural left, EcoRI right
        XCTAssertEqual(sortedFrags[0].start, 0)
        XCTAssertEqual(sortedFrags[0].end, 501)
        XCTAssertEqual(sortedFrags[0].length, 501)
        XCTAssertEqual(sortedFrags[0].leftEnd.overhangType, .blunt)
        XCTAssertEqual(sortedFrags[0].rightEnd.overhangType, .fivePrime)
        XCTAssertEqual(sortedFrags[0].rightEnd.overhangLen, 4)
        XCTAssertEqual(sortedFrags[0].rightEnd.overhangSeq5to3, "AATT")
        
        // Fragment 2: EcoRI left, PstI right
        XCTAssertEqual(sortedFrags[1].start, 501)
        XCTAssertEqual(sortedFrags[1].end, 1199)
        XCTAssertEqual(sortedFrags[1].length, 698)
        XCTAssertEqual(sortedFrags[1].leftEnd.overhangType, .fivePrime)
        XCTAssertEqual(sortedFrags[1].leftEnd.overhangLen, 4)
        XCTAssertEqual(sortedFrags[1].rightEnd.overhangType, .threePrime)
        XCTAssertEqual(sortedFrags[1].rightEnd.overhangLen, 4)
        
        // Fragment 3: PstI left, EcoRI right
        XCTAssertEqual(sortedFrags[2].start, 1199)
        XCTAssertEqual(sortedFrags[2].end, 1795)
        XCTAssertEqual(sortedFrags[2].length, 596)
        XCTAssertEqual(sortedFrags[2].leftEnd.overhangType, .threePrime)
        XCTAssertEqual(sortedFrags[2].leftEnd.overhangLen, 4)
        XCTAssertEqual(sortedFrags[2].rightEnd.overhangType, .fivePrime)
        XCTAssertEqual(sortedFrags[2].rightEnd.overhangLen, 4)
        
        // Fragment 4: EcoRI left, natural right
        XCTAssertEqual(sortedFrags[3].start, 1795)
        XCTAssertEqual(sortedFrags[3].end, 2000)
        XCTAssertEqual(sortedFrags[3].length, 205)
        XCTAssertEqual(sortedFrags[3].leftEnd.overhangType, .fivePrime)
        XCTAssertEqual(sortedFrags[3].leftEnd.overhangLen, 4)
        XCTAssertEqual(sortedFrags[3].rightEnd.overhangType, .blunt)
    }
    
    func test_incompatibleEnds_fail_compatibleConcrete() {
        // Test that incompatible pairs fail the compatibility check
        
        // Create different incompatible end types
        let ecoRIEnd = EndInfo(overhangType: .fivePrime,
                               overhangSeq5to3: "AATT",
                               sourceEnzyme: "EcoRI",
                               overhangLen: 4)
        
        let pstIEnd = EndInfo(overhangType: .threePrime,
                              overhangSeq5to3: "TGCA",
                              sourceEnzyme: "PstI",
                              overhangLen: 4)
        
        let bluntEnd = EndInfo(overhangType: .blunt,
                               overhangSeq5to3: nil,
                               sourceEnzyme: nil,
                               overhangLen: 0)
        
        let differentStickyEnd = EndInfo(overhangType: .fivePrime,
                                         overhangSeq5to3: "GATC",
                                         sourceEnzyme: "BamHI",
                                         overhangLen: 4)
        
        // Different overhang types should be incompatible
        XCTAssertFalse(compatibleConcrete(a: ecoRIEnd, b: pstIEnd),
                      "5' and 3' overhangs should be incompatible")
        
        // Sticky vs blunt should be incompatible
        XCTAssertFalse(compatibleConcrete(a: ecoRIEnd, b: bluntEnd),
                      "Sticky and blunt ends should be incompatible")
        
        // Different sticky sequences (same type) should be incompatible
        XCTAssertFalse(compatibleConcrete(a: ecoRIEnd, b: differentStickyEnd),
                      "Different sticky sequences should be incompatible")
        
        // Blunt-blunt without includeBlunt flag
        XCTAssertFalse(compatibleConcrete(a: bluntEnd, b: bluntEnd, includeBlunt: false),
                      "Blunt ends should be incompatible without includeBlunt flag")
        
        // Blunt-blunt with includeBlunt flag should be compatible
        XCTAssertTrue(compatibleConcrete(a: bluntEnd, b: bluntEnd, includeBlunt: true),
                     "Blunt ends should be compatible with includeBlunt flag")
    }
    
    func test_theoreticalCompatibility_NheI_XbaI() {
        // NheI and XbaI both create CTAG 5' overhangs (compatible)
        // NheI: GCTAGC, cut at 1 → G^CTAGC → 5' CTAG overhang
        // XbaI: TCTAGA, cut at 1 → T^CTAGA → 5' CTAG overhang
        
        let nheI = Enzyme(name: "NheI", site: "GCTAGC",
                          cutIndexTop: 1, cutIndexBottom: 5,
                          overhangType: .fivePrime, notes: nil)
        let xbaI = Enzyme(name: "XbaI", site: "TCTAGA",
                          cutIndexTop: 1, cutIndexBottom: 5,
                          overhangType: .fivePrime, notes: nil)
        
        let nheIEnd = theoreticalEnd(enzyme: nheI)
        let xbaIEnd = theoreticalEnd(enzyme: xbaI)
        
        // Both should produce CTAG overhang
        XCTAssertEqual(nheIEnd.template, "CTAG", "NheI should produce CTAG overhang")
        XCTAssertEqual(xbaIEnd.template, "CTAG", "XbaI should produce CTAG overhang")
        XCTAssertEqual(nheIEnd.type, .fivePrime)
        XCTAssertEqual(xbaIEnd.type, .fivePrime)
        XCTAssertEqual(nheIEnd.k, 4)
        XCTAssertEqual(xbaIEnd.k, 4)
        
        // They should be theoretically compatible
        XCTAssertTrue(compatibleTheoretical(a: nheIEnd, b: xbaIEnd),
                     "NheI and XbaI should be theoretically compatible (both produce CTAG)")
        
        // Test with actual digest to verify concrete compatibility
        let seqNheI = "AAAGCTAGCTT"
        let seqXbaI = "AAATCTAGATT"
        
        let engineNheI = DigestEngine(sequence: seqNheI, enzymes: [nheI])
        let engineXbaI = DigestEngine(sequence: seqXbaI, enzymes: [xbaI])
        
        let fragsNheI = engineNheI.digest(options: .init(circular: false))
        let fragsXbaI = engineXbaI.digest(options: .init(circular: false))
        
        let nheIStickyEnds = fragsNheI.flatMap { [$0.leftEnd, $0.rightEnd] }.filter { $0.overhangLen > 0 }
        let xbaIStickyEnds = fragsXbaI.flatMap { [$0.leftEnd, $0.rightEnd] }.filter { $0.overhangLen > 0 }
        
        XCTAssertFalse(nheIStickyEnds.isEmpty, "NheI should produce sticky ends")
        XCTAssertFalse(xbaIStickyEnds.isEmpty, "XbaI should produce sticky ends")
        
        // Verify the actual overhang sequences
        for end in nheIStickyEnds {
            XCTAssertEqual(end.overhangSeq5to3, "CTAG", "NheI should produce CTAG overhang")
        }
        for end in xbaIStickyEnds {
            XCTAssertEqual(end.overhangSeq5to3, "CTAG", "XbaI should produce CTAG overhang")
        }
        
        // Concrete compatibility test: CTAG is a palindrome, so it should be compatible
        XCTAssertTrue(compatibleConcrete(a: nheIStickyEnds[0], b: xbaIStickyEnds[0]),
                     "NheI and XbaI sticky ends should be concretely compatible")
    }
    
    func test_theoreticalIncompatibility() {
        // Test that clearly incompatible enzymes fail theoretical compatibility
        
        let ecoRI = Enzyme(name: "EcoRI", site: "GAATTC",
                           cutIndexTop: 1, cutIndexBottom: 5,
                           overhangType: .fivePrime, notes: nil)
        let pstI = Enzyme(name: "PstI", site: "CTGCAG",
                          cutIndexTop: 5, cutIndexBottom: 1,
                          overhangType: .threePrime, notes: nil)
        let bamHI = Enzyme(name: "BamHI", site: "GGATCC",
                           cutIndexTop: 1, cutIndexBottom: 5,
                           overhangType: .fivePrime, notes: nil)
        
        let ecoRIEnd = theoreticalEnd(enzyme: ecoRI)
        let pstIEnd = theoreticalEnd(enzyme: pstI)
        let bamHIEnd = theoreticalEnd(enzyme: bamHI)
        
        // Different overhang types (5' vs 3') should be incompatible
        XCTAssertFalse(compatibleTheoretical(a: ecoRIEnd, b: pstIEnd),
                      "EcoRI (5') and PstI (3') should be theoretically incompatible")
        
        // Same type but different sequences should be incompatible
        // EcoRI: AATT, BamHI: GATC
        XCTAssertFalse(compatibleTheoretical(a: ecoRIEnd, b: bamHIEnd),
                      "EcoRI (AATT) and BamHI (GATC) should be theoretically incompatible")
    }
    
    func test_boundaryPositions_circular() {
        // Test boundary positions for circular digest
        let ecoRI = Enzyme(name: "EcoRI", site: "GAATTC",
                           cutIndexTop: 1, cutIndexBottom: 5,
                           overhangType: .fivePrime, notes: nil)
        
        // Create a circular sequence with two EcoRI sites for proper testing
        // Site 1 at position 10, Site 2 at position 60
        let seq = String(repeating: "A", count: 10) + "GAATTC" + String(repeating: "T", count: 44) + "GAATTC" + String(repeating: "C", count: 34) // 100 bp total
        let engine = DigestEngine(sequence: seq, enzymes: [ecoRI])
        let frags = engine.digest(options: .init(circular: true))
        
        // Circular with two cuts should produce two fragments
        XCTAssertEqual(frags.count, 2, "Circular sequence with two cuts produces two fragments")
        
        // Expected cuts at positions 10+1=11 and 60+1=61
        // Fragment 1: [11, 61) = 50 bp
        // Fragment 2: [61, 11) wrapping = 39 bp (61 to 100 = 39, then 0 to 11 = 11, but wait...)
        // Actually: [61, 11+100) in extended form, which is 50 bp modulo
        
        let sortedFrags = frags.sorted { $0.start < $1.start }
        XCTAssertEqual(sortedFrags[0].start, 11)
        XCTAssertEqual(sortedFrags[0].end, 61)
        XCTAssertEqual(sortedFrags[0].length, 50)
        XCTAssertEqual(sortedFrags[0].leftEnd.overhangType, .fivePrime)
        XCTAssertEqual(sortedFrags[0].rightEnd.overhangType, .fivePrime)
        
        XCTAssertEqual(sortedFrags[1].start, 61)
        XCTAssertEqual(sortedFrags[1].end, 11)
        XCTAssertEqual(sortedFrags[1].length, 50)
        XCTAssertEqual(sortedFrags[1].leftEnd.overhangType, .fivePrime)
        XCTAssertEqual(sortedFrags[1].rightEnd.overhangType, .fivePrime)
    }
}
