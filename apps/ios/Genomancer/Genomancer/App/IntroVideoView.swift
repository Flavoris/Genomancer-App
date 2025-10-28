import SwiftUI
import AVFoundation
import UIKit

struct IntroVideoView: UIViewRepresentable {
    final class PlayerView: UIView {
        override static var layerClass: AnyClass { AVPlayerLayer.self }
        var playerLayer: AVPlayerLayer { layer as! AVPlayerLayer }
    }

    let player: AVPlayer
    let onFinished: () -> Void
    var videoGravity: AVLayerVideoGravity = .resizeAspect
    var loop: Bool = false

    func makeUIView(context: Context) -> PlayerView {
        let view = PlayerView()
        view.isOpaque = false
        view.backgroundColor = .clear
        view.playerLayer.player = player
        view.playerLayer.videoGravity = videoGravity
        view.playerLayer.isOpaque = false
        view.playerLayer.backgroundColor = UIColor.clear.cgColor

        // End observer
        if let item = player.currentItem {
            NotificationCenter.default.addObserver(
                context.coordinator,
                selector: #selector(Coordinator.itemDidEnd(_:)),
                name: .AVPlayerItemDidPlayToEndTime,
                object: item
            )
        }
        player.play()
        return view
    }

    func updateUIView(_ uiView: PlayerView, context: Context) {
        uiView.isOpaque = false
        uiView.backgroundColor = .clear
        uiView.playerLayer.player = player
        uiView.playerLayer.videoGravity = videoGravity
        uiView.playerLayer.isOpaque = false
        uiView.playerLayer.backgroundColor = UIColor.clear.cgColor
    }

    func makeCoordinator() -> Coordinator {
        Coordinator(onFinished: onFinished, player: player, loop: loop)
    }

    final class Coordinator: NSObject {
        private let onFinished: () -> Void
        private weak var player: AVPlayer?
        private let loop: Bool

        init(onFinished: @escaping () -> Void, player: AVPlayer, loop: Bool) {
            self.onFinished = onFinished
            self.player = player
            self.loop = loop
        }

        @objc func itemDidEnd(_ note: Notification) {
            print("🎥 Video player: AVPlayerItemDidPlayToEndTime notification received")
            guard let player = player else { 
                print("⚠️ Video player reference is nil")
                return 
            }
            if loop {
                print("🔁 Looping video")
                player.seek(to: .zero)
                player.play()
            } else {
                print("✅ Video finished, calling onFinished callback")
                onFinished()
            }
        }

        deinit {
            NotificationCenter.default.removeObserver(self)
        }
    }
}
