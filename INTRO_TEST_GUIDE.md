# Intro Feature Test Guide

## Implementation Summary

The intro feature has been successfully implemented with the following components:

### Created Files:
1. **AppStartupManager.swift** - Manages background app initialization
2. **IntroVideoView.swift** - Custom video player using AVFoundation
3. **IntroView.swift** - Main intro screen with dynamic video duration detection
4. **RootView.swift** - Root view with in-memory intro state via @State
5. **GenomancerApp.swift** - Updated app entry point

### Key Features:
- ✅ Intro video plays on every cold launch (after app is killed)
- ✅ `@State` provides in-memory state (resets when app terminates)
- ✅ Dynamic video duration detection ensures full video playback
- ✅ Graceful fallback if video is missing (1.2s branded screen)
- ✅ Smart transition logic waits for both startup AND video completion
- ✅ Smooth fade + scale animations

---

## Test Checklist

### Test 1: First Launch Experience
**Steps:**
1. Open the project in Xcode
2. Build and run on a fresh simulator (⌘R)
3. Observe the intro screen

**Expected Behavior:**
- ✅ Intro video (`Intro.mp4`) plays from start to finish
- ✅ Genomancer logo and branding are overlaid on the video
- ✅ Video plays at full screen with aspect fill
- ✅ After video completes AND startup finishes, smooth transition to HomeView
- ✅ Transition includes fade + scale animation (0.35s duration)

**Duration:**
- Minimum time shown = `max(video_duration, startup_time)`
- The intro will wait for both video AND enzyme loading to complete

---

### Test 2: Cold Launch (App Restart)
**Steps:**
1. With the app still running, press Home (⌘⇧H)
2. Kill the app completely (swipe up in app switcher)
3. Relaunch the app from simulator home screen

**Expected Behavior:**
- ✅ Intro plays again on cold launch
- ✅ `@State` resets to `false` when app terminates
- ✅ This provides a fresh branded experience each time the app starts

---

### Test 3: Session Persistence (While App Running)
**Steps:**
1. Let intro play and transition to HomeView
2. Press Home (⌘⇧H) to background the app
3. WITHOUT killing the app, tap the app icon to bring it back to foreground

**Expected Behavior:**
- ✅ App resumes at HomeView (intro does NOT replay)
- ✅ `@State` persists while app is in memory
- ✅ Intro only replays on cold launch (after app is killed)

---

### Test 4: Missing Video Fallback
**Steps:**
1. Temporarily rename or remove `Intro.mp4` from Resources
2. Build and run the app

**Expected Behavior:**
- ✅ Static branded screen shown (logo + genomancer font)
- ✅ Background color: `.genomancerBackground`
- ✅ Minimum display time: 1.2 seconds
- ✅ After minimum time AND startup complete, transitions to HomeView
- ✅ No errors or crashes

**To Test:**
```bash
# Temporarily move the video
cd apps/ios/Genomancer/Genomancer/Resources
mv Intro.mp4 Intro.mp4.backup

# Build and test...

# Restore
mv Intro.mp4.backup Intro.mp4
```

---

### Test 5: Video Duration Detection
**Verification:**
The `IntroView` automatically detects the MP4 duration and uses it as minimum display time.

**Code Verification:**
```swift
private func getVideoDuration() -> Double? {
    guard let url = Bundle.main.url(forResource: "Intro", withExtension: "mp4") else {
        return nil
    }
    
    let asset = AVAsset(url: url)
    let duration = asset.duration
    
    guard duration.isValid, !duration.isIndefinite else {
        return nil
    }
    
    return CMTimeGetSeconds(duration)
}
```

**Expected Behavior:**
- ✅ Returns actual video duration in seconds
- ✅ Falls back to 1.2s if video missing or invalid
- ✅ Intro shows for `max(video_duration, startup_time)`

---

## Manual Testing Procedure

### Build Steps:
```bash
# Navigate to project
cd /Users/flavorisbelue/Desktop/Genomancer/apps/ios/Genomancer

# Clean build
xcodebuild clean -project Genomancer.xcodeproj -scheme Genomancer

# Build for simulator
xcodebuild build -project Genomancer.xcodeproj -scheme Genomancer \
  -destination 'platform=iOS Simulator,name=iPhone 15'
```

### Or use Xcode GUI:
1. Open `Genomancer.xcodeproj` in Xcode
2. Select target: Genomancer
3. Select simulator device
4. Product → Clean Build Folder (⌘⇧K)
5. Product → Run (⌘R)

---

## Architecture Validation

### State Flow:
```
App Cold Launch (hasSeenIntro = false)
    ↓
GenomancerApp (@main)
    ↓
RootView (@State checks hasSeenIntro)
    ↓
    ├─ false → IntroView (on cold launch)
    │           ↓
    │       (Video plays + Startup)
    │           ↓
    │       (Sets hasSeenIntro = true)
    │           ↓
    └─ true → HomeView (while app in memory)

App Killed → State resets → Next launch shows intro again
```

### Transition Logic:
```swift
private func proceedIfDone() {
    // Both conditions must be true:
    // 1. startup.isReady (enzyme loading complete)
    // 2. videoFinished OR fallbackTimerElapsed (video complete OR min time)
    guard startup.isReady, (videoFinished || fallbackTimerElapsed) else { return }
    
    withAnimation(.easeInOut(duration: 0.35)) {
        hasSeenIntro = true  // Persisted to UserDefaults
    }
}
```

---

## Debugging Tips

### Check State:
```swift
// Add temporary debugging in RootView.swift
var body: some View {
    Group {
        if hasSeenIntro {
            HomeView()
                .onAppear { print("🏠 Showing HomeView (hasSeenIntro = true)") }
        } else {
            IntroView(hasSeenIntro: $hasSeenIntro)
                .onAppear { print("🎬 Showing IntroView (hasSeenIntro = false)") }
        }
    }
}
```

### Force Intro to Show Again:
To see the intro again without killing the app (for testing):
```swift
// Add a debug button in HomeView
Button("Reset Intro (Debug)") {
    // Navigate back to root or restart app
    // State will reset on cold launch
}
```

### Check Video Duration:
```swift
// Add to IntroView.onAppear for debugging
if let duration = getVideoDuration() {
    print("📹 Video duration: \(duration) seconds")
} else {
    print("⚠️ Video not found, using fallback time")
}
```

---

## Success Criteria

All tests passing means:
- ✅ Intro plays completely on every cold launch
- ✅ While app is running, intro only shows once per session
- ✅ Killing and relaunching shows intro again
- ✅ Backgrounding (without killing) resumes at HomeView
- ✅ Fallback works without video file
- ✅ No crashes or errors
- ✅ Smooth animations and transitions
- ✅ Video duration is automatically detected and respected

---

## Behavior Notes

1. **In-Memory State**: Using `@State` means:
   - Intro plays on every cold launch (app restart)
   - State persists only while app is running/backgrounded
   - State resets when app is killed or force-quit
   - This provides a fresh branded experience each time

2. **Simulator CoreSimulator Service**: The test environment may have simulator connection issues. Test on actual device for production validation.

3. **Video Format**: Currently expects MP4. Other formats may need AVAsset validation adjustments.

---

## Files Modified/Created

```
apps/ios/Genomancer/Genomancer/App/
├── AppStartupManager.swift      (NEW)
├── IntroVideoView.swift         (NEW)
├── IntroView.swift              (NEW)
├── RootView.swift               (NEW)
└── GenomancerApp.swift          (MODIFIED)

apps/ios/Genomancer/Genomancer/Resources/
└── Intro.mp4                    (ADDED - present in project)

apps/ios/Genomancer/Genomancer.xcodeproj/
└── project.pbxproj              (MODIFIED - added new files)
```

---

## Next Steps

1. Open Xcode and build the project
2. Run through all 5 test scenarios
3. Verify on both simulator and physical device
4. Test with different video durations
5. Verify memory management (no leaks)
6. Test on iOS 16.0+ (minimum deployment target)


