# Privacy Policy for Genomancer

**Last Updated: November 18, 2025**

## Our Commitment to Your Privacy

Genomancer is built with privacy as a core principle. We believe your DNA sequences and research data are yours alone. This privacy policy explains how we handle (or rather, don't handle) your data.

## The Short Version

**Genomancer collects absolutely no data. Period.**

Your sequences never leave your device. We don't track you, analyze you, or collect anything about you or your work.

---

## Detailed Privacy Policy

### 1. Information We Do NOT Collect

Genomancer does **not** collect, store, transmit, or process any of the following:

#### Personal Information
- ❌ No names, email addresses, or contact information
- ❌ No user accounts or authentication data
- ❌ No device identifiers or advertising IDs
- ❌ No location data or GPS coordinates
- ❌ No photos or camera access (except for document scanning if you use it)

#### Usage Data
- ❌ No analytics or usage statistics
- ❌ No crash reports or error logs
- ❌ No feature usage tracking
- ❌ No session duration or app interaction data
- ❌ No A/B testing or behavioral analysis

#### Research Data
- ❌ No DNA sequences you enter or import
- ❌ No enzyme selections or digest results
- ❌ No FASTA files or sequence data
- ❌ No fragment analysis or gel simulation data
- ❌ No plasmid maps or ligation results

#### Technical Data
- ❌ No IP addresses or network information
- ❌ No cookies or tracking pixels
- ❌ No diagnostic data or system logs
- ❌ No performance metrics

### 2. How Your Data is Processed

All data processing in Genomancer happens **100% locally on your device**:

#### Local-Only Processing
- ✅ DNA sequences are processed entirely in your iPhone/iPad's memory
- ✅ Enzyme matching algorithms run locally using the DigestCore Swift package
- ✅ Fragment calculations happen on-device with no cloud computation
- ✅ Gel simulations and plasmid maps are generated locally
- ✅ FASTA file parsing occurs entirely on your device

#### No Network Communication
- ✅ Genomancer does **not** require an internet connection to function
- ✅ The app makes **zero** network requests (except App Store checks by iOS)
- ✅ No data is uploaded to servers, cloud services, or third parties
- ✅ No remote APIs or external services are contacted

#### Data Storage
- ✅ Temporary data (current sequence, results) is stored in app memory only
- ✅ Data is cleared when you close the app or start a new analysis
- ✅ No persistent storage of sequences or results (unless you explicitly save)
- ✅ If you save sequences using iOS Share Sheet, that uses your device's standard file system (under iOS's control, not ours)

### 3. Third-Party Services

Genomancer uses **zero** third-party services:

- ❌ No analytics platforms (Google Analytics, Firebase, etc.)
- ❌ No advertising networks or ad providers
- ❌ No social media integrations
- ❌ No crash reporting services
- ❌ No customer support platforms
- ❌ No cloud storage or backup services
- ❌ No authentication providers

The only third-party code in Genomancer is:
- **Apple's iOS SDK** - Required for any iOS app (subject to [Apple's Privacy Policy](https://www.apple.com/legal/privacy/))
- **Swift Standard Library** - Part of the iOS platform

### 4. Permissions We Request

Genomancer may request the following iOS permissions, but only when needed:

#### File Access (when you import FASTA files)
- **Why:** To read DNA sequence files you explicitly choose to import
- **What we access:** Only the specific files you select via iOS file picker
- **What we do:** Read the file content into memory, parse the sequence, then discard the file
- **What we don't do:** We don't scan your file system, access other files, or store file metadata

#### Share Sheet Access (when you export results)
- **Why:** To let you save or share results using iOS's native sharing
- **What we access:** Whatever you choose to share (generated images, text, etc.)
- **What we do:** Provide data to iOS's Share Sheet for you to control
- **What we don't do:** We don't control where you send it or who receives it

### 5. Data Retention

**Genomancer retains no data.**

- App memory is cleared when you close the app
- iOS may keep app state for faster app switching (this is an iOS feature, not something we control)
- Killing the app completely clears all data from memory
- No data persists between app sessions unless you explicitly save files using iOS

### 6. Data Security

Since we don't collect or transmit data, there's no data to secure in transit or on servers. However:

- Your local sequences are protected by your device's security (passcode, Face ID, etc.)
- iOS's sandboxing prevents other apps from accessing Genomancer's memory
- All calculations happen in isolated app memory
- No credentials or API keys are stored (because there are none)

### 7. Children's Privacy

Genomancer is safe for users of all ages because:

- We collect no data from anyone, including children under 13
- The app has no user accounts, chat features, or social elements
- All content is educational and scientific
- We comply with COPPA (Children's Online Privacy Protection Act) by not collecting any personal information

### 8. International Users

Genomancer works the same way worldwide:

- No data collection means no cross-border data transfers
- No GDPR concerns (EU General Data Protection Regulation) because there's no personal data
- No CCPA concerns (California Consumer Privacy Act) because we don't sell or share data
- Compliant with all international privacy laws by design

### 9. Changes to This Privacy Policy

If we ever change how Genomancer handles data (which would be a major change):

- We will update this privacy policy with a new "Last Updated" date
- Significant changes will be noted in app update release notes
- We will never retroactively collect data from older app versions

You can always find the latest version of this policy at:
- On our GitHub repository

### 10. Your Rights

Under various privacy laws (GDPR, CCPA, etc.), you typically have rights to:

- **Access your data** → Not applicable (we have no data to access)
- **Delete your data** → Not applicable (we have no data to delete)
- **Port your data** → Your sequences stay on your device; you can export them anytime
- **Opt out of tracking** → Not applicable (we don't track)
- **Opt out of data sales** → Not applicable (we don't collect or sell data)

### 11. Contact Information

If you have questions about this privacy policy or Genomancer's privacy practices:

**Email:** flavorislab@gmail.com

**Response Time:** We typically respond within 3-5 business days

**GitHub Issues:** You can also open an issue on our GitHub repository

---

## Legal Compliance

### GDPR Compliance (European Union)
- **Lawful Basis:** Not applicable (no personal data collected)
- **Data Controller:** Not applicable (no data processing)
- **Data Processor:** Not applicable (no data processing)
- **DPO Contact:** Not required (no personal data)

### CCPA Compliance (California)
- **Data Sale:** We do not sell personal information
- **Data Sharing:** We do not share personal information
- **Data Collection:** We do not collect personal information
- **Opt-Out Rights:** Not applicable (nothing to opt out of)

### COPPA Compliance (Children's Privacy)
- **Parental Consent:** Not required (no personal information collected from children)
- **Verifiable Parental Consent:** Not applicable

### HIPAA (Health Information)
- **Note:** While DNA sequences might be considered health information in some contexts, Genomancer is designed for research and education, not medical diagnosis or treatment
- **PHI Handling:** We don't collect, store, or transmit any Protected Health Information
- **Business Associate Agreement:** Not applicable (no data transmission)

---

## Technical Privacy Details

### App Permissions (iOS Info.plist)
Genomancer requests minimal permissions:

```xml
<!-- File access for FASTA import -->
<key>NSDocumentDirectory</key>
<string>Import DNA sequence files (FASTA format)</string>

<!-- Photo library (if you save gel images) -->
<key>NSPhotoLibraryAddUsageDescription</key>
<string>Save gel simulation images to your photo library</string>
```

### Network Activity
- **Outbound connections:** None (0 requests)
- **Inbound connections:** None
- **Background network activity:** None
- **App Transport Security:** Not applicable (no network activity)

### Data Flow Diagram

```
User's Device
│
├─ Genomancer App (Sandboxed)
│  │
│  ├─ User enters DNA sequence → Stored in RAM only
│  ├─ DigestCore processes locally → Results in RAM only
│  ├─ UI displays results → Rendered locally
│  └─ User exports (optional) → iOS Share Sheet (user controls destination)
│
└─ No external servers, APIs, or cloud services
```

---

## Privacy By Design

Genomancer was built from the ground up with privacy as the foundation:

1. **Data Minimization:** We don't collect data because we don't need to
2. **Local Processing:** Everything happens on your device
3. **No Accounts:** No user accounts mean no authentication data
4. **No Analytics:** We trust users to provide feedback directly
5. **Open Source:** Transparency through code availability
6. **Offline First:** No internet required means no data leakage
7. **No Third Parties:** No SDKs that might collect data without our knowledge

---

## Frequently Asked Questions

### Q: Do you collect anonymous or aggregated usage data?
**A:** No. We collect nothing, not even anonymized data.

### Q: Do you use cookies?
**A:** No. Genomancer is a native iOS app, not a web app. There are no cookies.

### Q: What happens to my sequences when I close the app?
**A:** They're cleared from memory. Gone forever (unless you explicitly saved them to your device).

### Q: Can you recover my sequences if I lose them?
**A:** No. Since we never see or store your sequences, we can't recover them. Keep backups!

### Q: Does Apple collect data about Genomancer users?
**A:** Apple may collect aggregate data about App Store downloads and crashes through iOS's built-in systems. This is controlled by Apple, not by us. You can control this in iOS Settings → Privacy & Security.

### Q: Will you add analytics in the future?
**A:** We have no current plans to add analytics. If we ever did, it would be:
   - Clearly disclosed in an app update
   - Opt-in only
   - Privacy-preserving (anonymous, aggregated)
   - Subject to an updated privacy policy

### Q: How do I know you're telling the truth?
**A:** 
   - Use network monitoring tools to verify zero network activity
   - Check iOS Privacy settings → Genomancer (shows no permissions used)
   - Trust but verify - we encourage technical users to audit our claims

### Q: What about device fingerprinting?
**A:** We don't fingerprint devices. We don't even check your device model or iOS version (iOS provides this automatically to developers, but we don't use it).

---

## Conclusion

**Genomancer is built for scientists who value privacy.**

Your research is yours. Your sequences are yours. Your data stays with you.

If you have questions, concerns, or feedback about our privacy practices, please don't hesitate to contact us at **flavorislab@gmail.com**.

---

**Made with 🧬 and 🔒 for privacy-conscious researchers**

