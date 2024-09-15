import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

public class Bioinformatics {

    public static void main(String[] args) {
        int numSequences = 1000000; 
        int sequenceLength = 50; 
        double[] probabilities = {0.25, 0.25, 0.25, 0.25}; 
        int motifSize = 5; 
        double entropyThreshold = 1.5; 

        List<String> sequences = new ArrayList<>();
        SequenceGenerator generator = new SequenceGenerator(numSequences, sequenceLength, probabilities, entropyThreshold, sequences);
        generator.start();
        try {
            generator.join();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        saveSequencesToFile(sequences, "artificial_database.txt");

        MotifFinder finder = new MotifFinder(sequences, motifSize);
        finder.start();
        try {
            finder.join();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        String motif = finder.getBestMotif();
        System.out.println("The most frequent motif of size " + motifSize + " is: " + motif);
    }

    public static void saveSequencesToFile(List<String> sequences, String filename) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filename))) {
            for (String sequence : sequences) {
                writer.write(sequence);
                writer.newLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static double calculateEntropy(String sequence) {
        int[] counts = new int[4]; // A, C, G, T
        for (char base : sequence.toCharArray()) {
            switch (base) {
                case 'A':
                    counts[0]++;
                    break;
                case 'C':
                    counts[1]++;
                    break;
                case 'G':
                    counts[2]++;
                    break;
                case 'T':
                    counts[3]++;
                    break;
            }
        }

        double entropy = 0.0;
        int length = sequence.length();
        for (int count : counts) {
            if (count > 0) {
                double frequency = (double) count / length;
                entropy -= frequency * (Math.log(frequency) / Math.log(2));
            }
        }
        return entropy;
    }

    public static int countConsecutiveRepeats(String sequence) {
        int maxRepeats = 1;
        int currentRepeats = 1;
        for (int i = 1; i < sequence.length(); i++) {
            if (sequence.charAt(i) == sequence.charAt(i - 1)) {
                currentRepeats++;
                maxRepeats = Math.max(maxRepeats, currentRepeats);
            } else {
                currentRepeats = 1;
            }
        }
        return maxRepeats;
    }

    static class SequenceGenerator extends Thread {
        private int numSequences;
        private int sequenceLength;
        private double[] probabilities;
        private double entropyThreshold;
        private List<String> sequences;

        public SequenceGenerator(int numSequences, int sequenceLength, double[] probabilities, double entropyThreshold, List<String> sequences) {
            this.numSequences = numSequences;
            this.sequenceLength = sequenceLength;
            this.probabilities = probabilities;
            this.entropyThreshold = entropyThreshold;
            this.sequences = sequences;
        }

        @Override
        public void run() {
            for (int i = 0; i < numSequences; i++) {
                String sequence = generateSequence(sequenceLength, probabilities);
                if (calculateEntropy(sequence) >= entropyThreshold) {
                    synchronized (sequences) {
                        sequences.add(sequence);
                    }
                }
            }
        }

        private String generateSequence(int sequenceLength, double[] probabilities) {
            char[] nucleotides = {'A', 'C', 'G', 'T'};
            StringBuilder sequence = new StringBuilder(sequenceLength);
            Random random = new Random();

            for (int i = 0; i < sequenceLength; i++) {
                double rand = random.nextDouble();
                double cumulativeProbability = 0.0;
                for (int j = 0; j < probabilities.length; j++) {
                    cumulativeProbability += probabilities[j];
                    if (rand <= cumulativeProbability) {
                        sequence.append(nucleotides[j]);
                        break;
                    }
                }
            }
            return sequence.toString();
        }
    }

    static class MotifFinder extends Thread {
        private List<String> sequences;
        private int motifSize;
        private String bestMotif;
        private int maxCount;

        public MotifFinder(List<String> sequences, int motifSize) {
            this.sequences = sequences;
            this.motifSize = motifSize;
            this.bestMotif = "";
            this.maxCount = 0;
        }

        @Override
        public void run() {
            Map<String, Integer> motifCounts = new HashMap<>();
            for (String sequence : sequences) {
                for (int i = 0; i <= sequence.length() - motifSize; i++) {
                    String motif = sequence.substring(i, i + motifSize);
                    synchronized (motifCounts) {
                        motifCounts.put(motif, motifCounts.getOrDefault(motif, 0) + 1);
                    }
                }
            }

            for (Map.Entry<String, Integer> entry : motifCounts.entrySet()) {
                String motif = entry.getKey();
                int count = entry.getValue();
                if (count > maxCount || (count == maxCount && countConsecutiveRepeats(motif) > countConsecutiveRepeats(bestMotif))) {
                    maxCount = count;
                    bestMotif = motif;
                }
            }
        }

        public String getBestMotif() {
            return bestMotif;
        }
    }
}