import SwiftUI
import AVFoundation

<<<<<<< ours
<<<<<<< ours
struct IntroVideoView: UIViewRepresentable {
=======
struct TransparentVideoView: UIViewRepresentable {
>>>>>>> theirs
=======
struct IntroVideoView: UIViewRepresentable {
>>>>>>> theirs
    final class PlayerView: UIView {
        override static var layerClass: AnyClass { AVPlayerLayer.self }
        var playerLayer: AVPlayerLayer { layer as! AVPlayerLayer }
    }

    let player: AVPlayer
    let onFinished: () -> Void
    var videoGravity: AVLayerVideoGravity = .resizeAspectFill
    var loop: Bool = false

    func makeUIView(context: Context) -> PlayerView {
        let view = PlayerView()
        view.isOpaque = false
        view.backgroundColor = .clear
<<<<<<< ours
<<<<<<< ours
        view.playerLayer.player = player
        view.playerLayer.videoGravity = videoGravity
        view.playerLayer.backgroundColor = UIColor(Color.genomancerBackground).cgColor

        // End observer
=======

        let layer = view.playerLayer
        layer.isOpaque = false
        layer.backgroundColor = UIColor.clear.cgColor
        layer.videoGravity = videoGravity
        layer.player = player
        layer.frame = UIScreen.main.bounds

>>>>>>> theirs
=======
        view.playerLayer.player = player
        view.playerLayer.videoGravity = videoGravity
        view.playerLayer.backgroundColor = UIColor(Color.genomancerBackground).cgColor

        // End observer
>>>>>>> theirs
        if let item = player.currentItem {
            NotificationCenter.default.addObserver(
                context.coordinator,
                selector: #selector(Coordinator.itemDidEnd(_:)),
                name: .AVPlayerItemDidPlayToEndTime,
                object: item
            )
        }
<<<<<<< ours
<<<<<<< ours
=======

        player.isMuted = true
>>>>>>> theirs
=======
>>>>>>> theirs
        player.play()
        return view
    }

    func updateUIView(_ uiView: PlayerView, context: Context) {
<<<<<<< ours
<<<<<<< ours
        uiView.playerLayer.player = player
        uiView.playerLayer.videoGravity = videoGravity
        uiView.playerLayer.backgroundColor = UIColor(Color.genomancerBackground).cgColor
=======
        let layer = uiView.playerLayer
        layer.player = player
        layer.videoGravity = videoGravity
        layer.frame = UIScreen.main.bounds
>>>>>>> theirs
=======
        uiView.playerLayer.player = player
        uiView.playerLayer.videoGravity = videoGravity
        uiView.playerLayer.backgroundColor = UIColor(Color.genomancerBackground).cgColor
>>>>>>> theirs
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
<<<<<<< ours
<<<<<<< ours
            guard let player = player else { 
                print("⚠️ Video player reference is nil")
                return 
            }
=======
            guard let player = player else {
=======
            guard let player = player else { 
>>>>>>> theirs
                print("⚠️ Video player reference is nil")
                return 
            }
<<<<<<< ours

>>>>>>> theirs
=======
>>>>>>> theirs
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
<<<<<<< ours
<<<<<<< ours

=======
>>>>>>> theirs
=======

>>>>>>> theirs
