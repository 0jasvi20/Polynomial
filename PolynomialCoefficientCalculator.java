import java.math.BigInteger;
import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class PolynomialCoefficientCalculator {
    
    static class Point {
        int x;
        BigInteger y;
        
        Point(int x, BigInteger y) {
            this.x = x;
            this.y = y;
        }
        
        @Override
        public String toString() {
            return "(" + x + ", " + y + ")";
        }
    }
    
    static class TestResult {
        String testName;
        BigInteger constantC;
        boolean success;
        String error;
        
        TestResult(String testName, BigInteger constantC) {
            this.testName = testName;
            this.constantC = constantC;
            this.success = true;
        }
        
        TestResult(String testName, String error) {
            this.testName = testName;
            this.error = error;
            this.success = false;
        }
    }
    
    public static void main(String[] args) {
        List<TestResult> results = new ArrayList<>();
        
        // Test Case 1 - First sample test case
        String testCase1 = """
            {
              "keys": {
                "n": 4,
                "k": 3
              },
              "1": {
                "base": "10",
                "value": "4"
              },
              "2": {
                "base": "2",
                "value": "111"
              },
              "3": {
                "base": "10",
                "value": "12"
              },
              "6": {
                "base": "4",
                "value": "213"
              }
            }
            """;
        
        // Test Case 2 - Second sample test case (large numbers)
        String testCase2 = """
            {
              "keys": {
                "n": 10,
                "k": 7
              },
              "1": {
                "base": "6",
                "value": "13444211440455345511"
              },
              "2": {
                "base": "15",
                "value": "aed7015a346d635"
              },
              "3": {
                "base": "15",
                "value": "6aeeb69631c227c"
              },
              "4": {
                "base": "16",
                "value": "e1b5e05623d881f"
              },
              "5": {
                "base": "8",
                "value": "316034514573652620673"
              },
              "6": {
                "base": "3",
                "value": "2122212201122002221120200210011020220200"
              },
              "7": {
                "base": "3",
                "value": "20120221122211000100210021102001201112121"
              },
              "8": {
                "base": "6",
                "value": "20220554335330240002224253"
              },
              "9": {
                "base": "12",
                "value": "45153788322a1255483"
              },
              "10": {
                "base": "7",
                "value": "1101613130313526312514143"
              }
            }
            """;
        
        // Process both test cases
        System.out.println("PROCESSING ALL TEST CASES...\n");
        
        results.add(processTestCase("Sample Test Case 1", testCase1, true));
        results.add(processTestCase("Sample Test Case 2 (Large)", testCase2, true));
        
        // Print final summary
        printSummary(results);
    }
    
    public static TestResult processTestCase(String testName, String jsonString, boolean showDetails) {
        System.out.println("\n" + "=".repeat(70));
        System.out.println("PROCESSING " + testName.toUpperCase());
        System.out.println("=".repeat(70));
        
        try {
            BigInteger constantC = findConstantCoefficient(jsonString, showDetails);
            System.out.println("\n✅ SUCCESS - Constant coefficient 'c': " + constantC);
            return new TestResult(testName, constantC);
        } catch (Exception e) {
            System.err.println("\n❌ ERROR - " + e.getMessage());
            e.printStackTrace();
            return new TestResult(testName, e.getMessage());
        }
    }
    
    public static void printSummary(List<TestResult> results) {
        System.out.println("\n" + "=".repeat(70));
        System.out.println("FINAL SUMMARY - ALL TEST CASES");
        System.out.println("=".repeat(70));
        
        for (TestResult result : results) {
            if (result.success) {
                System.out.println("✅ " + result.testName + ":");
                System.out.println("   Constant coefficient 'c' = " + result.constantC);
            } else {
                System.out.println("❌ " + result.testName + ": FAILED");
                System.out.println("   Error: " + result.error);
            }
            System.out.println();
        }
        
        long successCount = results.stream().mapToLong(r -> r.success ? 1 : 0).sum();
        System.out.println("RESULTS: " + successCount + "/" + results.size() + " test cases passed");
    }
    
    public static BigInteger findConstantCoefficient(String jsonString, boolean showDetails) {
        // Parse JSON manually
        int n = parseIntValue(jsonString, "\"n\":\\s*(\\d+)");
        int k = parseIntValue(jsonString, "\"k\":\\s*(\\d+)");
        
        if (showDetails) {
            System.out.println("Parameters:");
            System.out.println("  n (total roots provided): " + n);
            System.out.println("  k (minimum roots needed): " + k);
            System.out.println("  Polynomial degree: " + (k-1));
            System.out.println();
        }
        
        // Extract and convert points - search for ALL numeric keys, not just 1 to n
        List<Point> points = new ArrayList<>();
        
        // Find all point patterns in the JSON
        Pattern pointPattern = Pattern.compile("\"(\\d+)\":\\s*\\{[^}]*\"base\":\\s*\"(\\d+)\"[^}]*\"value\":\\s*\"([^\"]+)\"[^}]*\\}");
        Matcher matcher = pointPattern.matcher(jsonString);
        
        while (matcher.find()) {
            try {
                int x = Integer.parseInt(matcher.group(1));
                int base = Integer.parseInt(matcher.group(2));
                String valueStr = matcher.group(3);
                
                // Convert from given base to decimal using BigInteger
                BigInteger y = new BigInteger(valueStr, base);
                
                points.add(new Point(x, y));
                if (showDetails) {
                    System.out.println("Point " + x + ": " + valueStr + " (base " + base + ")");
                    System.out.println("         = " + y + " (decimal)");
                }
            } catch (NumberFormatException e) {
                System.err.println("Error parsing point: " + e.getMessage());
            }
        }
        
        if (showDetails) {
            System.out.println("\nFound " + points.size() + " points total");
        }
        
        // Sort points by x value for consistent processing
        points.sort((p1, p2) -> Integer.compare(p1.x, p2.x));
        
        // Take first k points for calculation
        List<Point> selectedPoints = points.subList(0, Math.min(k, points.size()));
        if (showDetails) {
            System.out.println("Using first " + selectedPoints.size() + " points for calculation:");
            selectedPoints.forEach(p -> System.out.println("  " + p));
            System.out.println();
        }
        
        // Calculate using Lagrange interpolation
        if (showDetails) {
            System.out.println("=== LAGRANGE INTERPOLATION ===");
        }
        BigInteger result = lagrangeInterpolationAtZero(selectedPoints, showDetails);
        
        // Demonstrate Sridharacharya formula on first 3 points if possible
        if (showDetails && selectedPoints.size() >= 3) {
            System.out.println("\n=== SRIDHARACHARYA FORMULA DEMONSTRATION ===");
            System.out.println("(Applied to quadratic formed by first 3 points)");
            demonstrateQuadraticReconstruction(selectedPoints.subList(0, 3));
        }
        
        return result;
    }
    
    /**
     * Parse integer value from JSON string using regex
     */
    private static int parseIntValue(String json, String pattern) {
        Pattern p = Pattern.compile(pattern);
        Matcher m = p.matcher(json);
        if (m.find()) {
            return Integer.parseInt(m.group(1));
        }
        throw new RuntimeException("Could not find pattern: " + pattern);
    }
    
    /**
     * Demonstrate quadratic reconstruction and Sridharacharya formula
     */
    private static void demonstrateQuadraticReconstruction(List<Point> points) {
        System.out.println("Reconstructing quadratic from points: " + points);
        
        try {
            // Use BigDecimal for matrix operations
            BigDecimal[][] matrix = new BigDecimal[3][4];
            MathContext mc = new MathContext(50, RoundingMode.HALF_UP);
            
            for (int i = 0; i < 3; i++) {
                Point p = points.get(i);
                matrix[i][0] = new BigDecimal(p.x * p.x); // x²
                matrix[i][1] = new BigDecimal(p.x);        // x
                matrix[i][2] = BigDecimal.ONE;             // constant term
                matrix[i][3] = new BigDecimal(p.y);        // y value
            }
            
            BigDecimal[] coeffs = gaussianEliminationBigDecimal(matrix, mc);
            double a = coeffs[0].doubleValue();
            double b = coeffs[1].doubleValue();  
            double c = coeffs[2].doubleValue();
            
            System.out.printf("Reconstructed quadratic: %.6fx² + %.6fx + %.6f%n", a, b, c);
            
            // Apply Sridharacharya formula
            applySridharacharyaFormula(a, b, c);
        } catch (Exception e) {
            System.err.println("Error in quadratic reconstruction: " + e.getMessage());
        }
    }
    
    /**
     * Apply Sridharacharya Formula (Quadratic Formula) to find roots
     * Formula: x = (-b ± √(b² - 4ac)) / (2a)
     */
    private static void applySridharacharyaFormula(double a, double b, double c) {
        System.out.println("\n--- SRIDHARACHARYA FORMULA APPLICATION ---");
        System.out.printf("Quadratic equation: %.6fx² + %.6fx + %.6f = 0%n", a, b, c);
        
        if (Math.abs(a) < 1e-15) {
            System.out.println("Not a quadratic (a ≈ 0)");
            if (Math.abs(b) > 1e-15) {
                double root = -c / b;
                System.out.printf("Linear root: x = %.6f%n", root);
            }
            return;
        }
        
        // Calculate discriminant: b² - 4ac
        double discriminant = b * b - 4 * a * c;
        System.out.printf("Discriminant = b² - 4ac = (%.6f)² - 4(%.6f)(%.6f) = %.6f%n", b, a, c, discriminant);
        
        if (discriminant > 0) {
            double sqrt_d = Math.sqrt(discriminant);
            double root1 = (-b + sqrt_d) / (2 * a);
            double root2 = (-b - sqrt_d) / (2 * a);
            
            System.out.println("Two distinct real roots:");
            System.out.printf("x₁ = (-%.6f + √%.6f) / %.6f = %.6f%n", b, discriminant, 2*a, root1);
            System.out.printf("x₂ = (-%.6f - √%.6f) / %.6f = %.6f%n", b, discriminant, 2*a, root2);
            
        } else if (Math.abs(discriminant) < 1e-15) {
            double root = -b / (2 * a);
            System.out.printf("One repeated real root: x = %.6f%n", root);
            
        } else {
            double realPart = -b / (2 * a);
            double imagPart = Math.sqrt(-discriminant) / (2 * a);
            System.out.printf("Two complex conjugate roots:%n");
            System.out.printf("x₁ = %.6f + %.6fi%n", realPart, imagPart);
            System.out.printf("x₂ = %.6f - %.6fi%n", realPart, imagPart);
        }
    }
    
    /**
     * Lagrange interpolation to find polynomial value at x=0 (constant term)
     * This is the main method to find coefficient 'c'
     */
    private static BigInteger lagrangeInterpolationAtZero(List<Point> points, boolean showDetails) {
        if (showDetails) {
            System.out.println("Computing polynomial value at x=0 using Lagrange interpolation:");
            System.out.println("f(0) = Σ yᵢ × Lᵢ(0), where Lᵢ(0) = ∏(0-xⱼ)/(xᵢ-xⱼ) for j≠i");
            System.out.println();
        }
        
        // Use exact rational arithmetic
        BigInteger resultNumerator = BigInteger.ZERO;
        BigInteger resultDenominator = BigInteger.ONE;
        
        for (int i = 0; i < points.size(); i++) {
            Point pi = points.get(i);
            
            // Calculate Lagrange basis polynomial Lᵢ(0)
            BigInteger numerator = BigInteger.ONE;
            BigInteger denominator = BigInteger.ONE;
            
            for (int j = 0; j < points.size(); j++) {
                if (i != j) {
                    Point pj = points.get(j);
                    // Lᵢ(0) term: (0 - xⱼ) / (xᵢ - xⱼ) = -xⱼ / (xᵢ - xⱼ)
                    numerator = numerator.multiply(BigInteger.valueOf(-pj.x));
                    denominator = denominator.multiply(BigInteger.valueOf(pi.x - pj.x));
                }
            }
            
            if (showDetails) {
                System.out.printf("L%d(0) = %s/%s%n", i, numerator, denominator);
            }
            
            // Add yᵢ × Lᵢ(0) to result using rational arithmetic
            BigInteger termNumerator = pi.y.multiply(numerator).multiply(resultDenominator);
            BigInteger newDenominator = denominator.multiply(resultDenominator);
            
            resultNumerator = resultNumerator.multiply(denominator).add(termNumerator);
            resultDenominator = newDenominator;
            
            // Reduce fraction by GCD
            BigInteger gcd = resultNumerator.gcd(resultDenominator);
            resultNumerator = resultNumerator.divide(gcd);
            resultDenominator = resultDenominator.divide(gcd);
            
            if (showDetails) {
                System.out.printf("After adding term %d: %s/%s%n", i, resultNumerator, resultDenominator);
            }
        }
        
        // Final result should be an integer for polynomial coefficients
        if (!resultDenominator.equals(BigInteger.ONE)) {
            if (showDetails) {
                System.out.printf("Warning: Non-integer result %s/%s, performing integer division%n", 
                                resultNumerator, resultDenominator);
            }
            return resultNumerator.divide(resultDenominator);
        }
        
        return resultNumerator;
    }
    
    /**
     * Gaussian elimination using BigDecimal for better precision
     */
    private static BigDecimal[] gaussianEliminationBigDecimal(BigDecimal[][] matrix, MathContext mc) {
        int n = matrix.length;
        
        // Forward elimination
        for (int i = 0; i < n; i++) {
            // Partial pivoting
            int maxRow = i;
            for (int k = i + 1; k < n; k++) {
                if (matrix[k][i].abs().compareTo(matrix[maxRow][i].abs()) > 0) {
                    maxRow = k;
                }
            }
            
            // Swap rows
            BigDecimal[] temp = matrix[i];
            matrix[i] = matrix[maxRow];
            matrix[maxRow] = temp;
            
            // Eliminate column
            for (int k = i + 1; k < n; k++) {
                if (matrix[i][i].compareTo(BigDecimal.ZERO) != 0) {
                    BigDecimal factor = matrix[k][i].divide(matrix[i][i], mc);
                    for (int j = i; j <= n; j++) {
                        matrix[k][j] = matrix[k][j].subtract(factor.multiply(matrix[i][j], mc), mc);
                    }
                }
            }
        }
        
        // Back substitution
        BigDecimal[] solution = new BigDecimal[n];
        for (int i = n - 1; i >= 0; i--) {
            solution[i] = matrix[i][n];
            for (int j = i + 1; j < n; j++) {
                solution[i] = solution[i].subtract(matrix[i][j].multiply(solution[j], mc), mc);
            }
            if (matrix[i][i].compareTo(BigDecimal.ZERO) != 0) {
                solution[i] = solution[i].divide(matrix[i][i], mc);
            }
        }
        
        return solution;
    }
}