import SwiftUI
import AVFoundation

struct IntroView: View {
    @Binding var hasSeenIntro: Bool
    @StateObject private var startup = AppStartupManager()

    @State private var videoFinished = false

    var body: some View {
        ZStack {
            // Background color
            Color.genomancerBackground
                .ignoresSafeArea()
            
            // Centered video
            if let url = Bundle.main.url(forResource: "Intro", withExtension: "mov") {
                VideoBackground(url: url) {
                    print("🎬 Video finished callback received")
                    videoFinished = true
                    proceedIfDone()
                }
                .onAppear {
                    print("📹 Video URL found: \(url.lastPathComponent)")
                }
            } else {
                Text("Video not found")
                    .foregroundColor(.white)
                    .onAppear {
                        print("❌ Intro.mov not found in bundle")
                    }
            }
        }
        .onAppear {
            print("🚀 IntroView appeared")
            startup.performStartup()
        }
        .onChange(of: startup.isReady) { _ in proceedIfDone() }
        .onChange(of: videoFinished) { _ in proceedIfDone() }
    }

    private func proceedIfDone() {
        // Proceed ONLY when both startup finished AND video finished
        // This ensures the entire video plays before transitioning
        guard startup.isReady, videoFinished else { 
            // Debug logging
            if !startup.isReady {
                print("⏳ Waiting for startup to complete...")
            }
            if !videoFinished {
                print("⏳ Waiting for video to finish playing...")
            }
            return 
        }
        
        print("✅ Intro complete! Startup ready: \(startup.isReady), Video finished: \(videoFinished)")
        print("✅ Transitioning to HomeView")
        withAnimation(.easeInOut(duration: 0.35)) {
            hasSeenIntro = true
        }
    }
}

private struct VideoBackground: View {
    let url: URL
    let onFinished: () -> Void
    
    @State private var player: AVPlayer?

    var body: some View {
        Group {
            if let player = player {
                IntroVideoView(player: player, onFinished: onFinished, videoGravity: .resizeAspect, loop: false)
                    .frame(maxWidth: .infinity, maxHeight: .infinity)
            } else {
                Text("Loading video...")
                    .foregroundColor(.white)
            }
        }
        .onAppear {
            print("🎥 VideoBackground appeared, creating player")
            let newPlayer = AVPlayer(url: url)
            newPlayer.volume = 1.0
            print("🎥 Player created, starting playback")
            self.player = newPlayer
            newPlayer.play()
        }
    }
}
