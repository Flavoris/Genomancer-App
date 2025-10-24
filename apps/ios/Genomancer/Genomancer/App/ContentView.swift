import SwiftUI
import DigestCore

struct ContentView: View {
    @State private var digestCore = DigestCore()
    @State private var sequence = "ATCGATCGATCG"
    @State private var analysisResult = ""
    
    var body: some View {
        NavigationView {
            VStack(spacing: 20) {
                Text("Genomancer")
                    .font(.largeTitle)
                    .fontWeight(.bold)
                    .foregroundColor(.genomancerText)
                    .padding()
                
                Text("DNA Digestion Analysis Tool")
                    .font(.headline)
                    .foregroundColor(.genomancerText)
                
                VStack(alignment: .leading, spacing: 10) {
                    Text("DNA Sequence:")
                        .font(.headline)
                        .foregroundColor(.genomancerText)
                    
                    TextField("Enter DNA sequence", text: $sequence)
                        .textFieldStyle(RoundedBorderTextFieldStyle())
                        .autocapitalization(.allCharacters)
                        .disableAutocorrection(true)
                }
                .padding()
                
                Button("Analyze Sequence") {
                    analysisResult = digestCore.analyzeSequence(sequence)
                }
                .buttonStyle(.borderedProminent)
                .disabled(sequence.isEmpty)
                
                if !analysisResult.isEmpty {
                    VStack(alignment: .leading, spacing: 10) {
                        Text("Analysis Result:")
                            .font(.headline)
                            .foregroundColor(.genomancerText)
                        
                        Text(analysisResult)
                            .padding()
                            .background(Color.genomancerCardBackground)
                            .cornerRadius(8)
                    }
                    .padding()
                }
                
                Spacer()
            }
            .background(Color.genomancerBackground)
            .navigationTitle("Genomancer")
            .navigationBarTitleDisplayMode(.inline)
        }
    }
}

struct ContentView_Previews: PreviewProvider {
    static var previews: some View {
        ContentView()
    }
}
