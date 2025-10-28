import SwiftUI

struct RootView: View {
    // In-memory state: resets on cold launch (after app is killed)
    @State private var hasSeenIntro: Bool = false

    var body: some View {
        Group {
            if hasSeenIntro {
                HomeView()
                    .transition(.opacity.combined(with: .scale))
            } else {
                IntroView(hasSeenIntro: $hasSeenIntro)
                    .transition(.opacity)
            }
        }
    }
}

