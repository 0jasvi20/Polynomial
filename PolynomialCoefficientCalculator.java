import org.json.JSONObject;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

public class PolynomialFromJson {

    // Method to convert value to decimal based on base
    public static int convertToDecimal(String value, int base) {
        return Integer.parseInt(value, base);
    }

    // Method to calculate polynomial coefficients for quadratic
    public static void calculateQuadraticCoefficients(double r1, double r2) {
        double a = 1.0;
        double b = -(r1 + r2);
        double c = r1 * r2;

        System.out.println("\n--- Polynomial Coefficients ---");
        System.out.println("a = " + a);
        System.out.println("b = " + b);
        System.out.println("c = " + c);

        // Verify roots using Sridharacharya formula
        double discriminant = b * b - 4 * a * c;
        if (discriminant < 0) {
            System.out.println("The polynomial has complex roots.");
        } else {
            double root1 = (-b + Math.sqrt(discriminant)) / (2 * a);
            double root2 = (-b - Math.sqrt(discriminant)) / (2 * a);
            System.out.println("\n--- Verification ---");
            System.out.println("Calculated Root 1: " + root1);
            System.out.println("Calculated Root 2: " + root2);
        }
    }

    public static void main(String[] args) {
        try {
            // Read JSON from file
            String content = new String(Files.readAllBytes(Paths.get("data.json")));
            JSONObject jsonObject = new JSONObject(content);

            // Extract n and k
            JSONObject keys = jsonObject.getJSONObject("keys");
            int n = keys.getInt("n");
            int k = keys.getInt("k");
            int degree = k - 1;

            System.out.println("Total roots provided (n): " + n);
            System.out.println("Minimum roots required (k): " + k);
            System.out.println("Degree of polynomial (m): " + degree);

            // Extract all roots into a list
            List<Integer> roots = new ArrayList<>();
            for (String key : jsonObject.keySet()) {
                if (!key.equals("keys")) {
                    JSONObject entry = jsonObject.getJSONObject(key);
                    int base = Integer.parseInt(entry.getString("base"));
                    String value = entry.getString("value");
                    roots.add(convertToDecimal(value, base));
                }
            }

            // Sort roots to have a deterministic order
            Collections.sort(roots);
            System.out.println("\nAll roots (decimal): " + roots);

            // Check if we have enough roots
            if (roots.size() < degree) {
                System.out.println("Not enough roots to form a polynomial of degree " + degree);
                return;
            }

            // For quadratic (degree=2), take first 2 roots
            if (degree == 2) {
                double r1 = roots.get(0);
                double r2 = roots.get(1);
                calculateQuadraticCoefficients(r1, r2);
            } else {
                System.out.println("Polynomial of degree " + degree + " is not implemented yet.");
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
