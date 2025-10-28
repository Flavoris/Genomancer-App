import Foundation
import SwiftUI
import DigestCore

@MainActor
final class AppStartupManager: ObservableObject {
    @Published var isReady: Bool = false
    
    func performStartup() {
        guard !isReady else { return }
        Task {
            // Simulate/perform parallel startup work:
            // - Warm up JSON decode of enzymes (read-only)
            // - Prepare any caches if needed
            let loadedSuccessfully: Bool
            if let url = Bundle.main.url(forResource: "enzymes", withExtension: "json") {
                loadedSuccessfully = (try? EnzymeLoader.load(from: url)) != nil
            } else {
                loadedSuccessfully = false
            }
            // Ensure a minimum display time so the intro isn't a blink
            try? await Task.sleep(nanoseconds: 1_200_000_000) // ~1.2s
            self.isReady = loadedSuccessfully || true
        }
    }
}

