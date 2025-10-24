import XCTest
@testable import DigestCore

final class GelLayoutTests: XCTestCase {
    func testTopAtMaxBP() {
        let top: CGFloat = 10
        let height: CGFloat = 200
        let minBP = 100
        let maxBP = 5000
        // Use reference percent (2.0%) for predictable positioning
        let y = GelLayout.y(bp: maxBP, minBP: minBP, maxBP: maxBP, top: top, height: height, gelPercent: 2.0)
        XCTAssertEqual(y, top, accuracy: 0.5)
    }

    func testBottomAtMinBP() {
        let top: CGFloat = 10
        let height: CGFloat = 200
        let minBP = 100
        let maxBP = 5000
        // Use reference percent (2.0%) for predictable positioning
        let y = GelLayout.y(bp: minBP, minBP: minBP, maxBP: maxBP, top: top, height: height, gelPercent: 2.0)
        XCTAssertEqual(y, top + height, accuracy: 0.5)
    }

    func testMonotonicity() {
        let top: CGFloat = 0
        let height: CGFloat = 300
        let minBP = 100
        let maxBP = 10000
        let values = [120, 200, 500, 1000, 2000, 5000, 9000]
        let ys = values.map { GelLayout.y(bp: $0, minBP: minBP, maxBP: maxBP, top: top, height: height, gelPercent: 1.5) }
        // Larger bp should map higher (smaller y)
        for i in 1..<ys.count {
            XCTAssertLessThanOrEqual(ys[i], ys[i-1] + 1e-6)
        }
    }

    func testLogLinearityInPercent() {
        let bp = 800
        let Cs = [1.0, 1.5, 2.0, 2.5, 3.0]
        let ys = Cs.map { log(GelLayout.mobility(bp: bp, gelPercent: $0)) }
        // successive differences ~ constant
        let diffs = zip(ys.dropFirst(), ys).map(-)
        let mean = diffs.reduce(0, +) / Double(diffs.count)
        let maxDev = diffs.map { abs($0 - mean) }.max() ?? 0
        XCTAssertLessThan(maxDev, 0.08)
    }

    func testNoTopClumpingInY() {
        let C = 2.0
        let top: CGFloat = 0, h: CGFloat = 1000
        let sizes = [1200, 900, 700, 500, 300, 150, 80]
        let ys = sizes.map { GelLayout.y(bp: $0, minBP: 80, maxBP: 1200, top: top, height: h, gelPercent: C) }
        // Strictly increasing with size decreasing
        for i in 1..<ys.count { XCTAssertLessThan(ys[i-1], ys[i]) }
        // Upper-half should not contain all large fragments (avoid clumping):
        let topRegion = top + h * 0.15
        let numInTop15 = ys.filter { $0 < topRegion }.count
        XCTAssertLessThan(numInTop15, 3, "Too many bands collapsed near the top")
    }

    func testMobilityMonotonicWithSize() {
        let C = 2.0
        let m100  = GelLayout.mobility(bp: 100,  gelPercent: C)
        let m300  = GelLayout.mobility(bp: 300,  gelPercent: C)
        let m1000 = GelLayout.mobility(bp: 1000, gelPercent: C)
        let m3000 = GelLayout.mobility(bp: 3000, gelPercent: C)
        XCTAssertTrue(m100  > m300)
        XCTAssertTrue(m300  > m1000)
        XCTAssertTrue(m1000 > m3000)
    }

    func testMobilityDecreasesWithGelPercent() {
        let bp = 1000
        let m1 = GelLayout.mobility(bp: bp, gelPercent: 1.0)
        let m2 = GelLayout.mobility(bp: bp, gelPercent: 2.0)
        let m3 = GelLayout.mobility(bp: bp, gelPercent: 3.0)
        XCTAssertTrue(m1 > m2 && m2 > m3)
    }

    func testYMappingDirection() {
        let top: CGFloat = 0
        let height: CGFloat = 300
        let ySmall = GelLayout.y(bp: 100,  minBP: 50,  maxBP: 5000, top: top, height: height, gelPercent: 2.0)
        let yLarge = GelLayout.y(bp: 3000, minBP: 50,  maxBP: 5000, top: top, height: height, gelPercent: 2.0)
        XCTAssertTrue(ySmall > yLarge, "Smaller fragments should appear lower (greater y).")
    }
}


