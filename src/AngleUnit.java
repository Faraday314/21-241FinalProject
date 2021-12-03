import org.jetbrains.annotations.NotNull;

import java.util.function.Function;

/**
 * A enum used for convntly
 */
public enum AngleUnit {
    DEGREES, RADIANS;

    public @NotNull Function<Double, Double> convertTo(AngleUnit angleUnit) {
        if(this.equals(angleUnit)) {
            return Double::doubleValue;
        }
        else if(angleUnit.equals(RADIANS)) {
            return Math::toRadians;
        }
        else {
            return Math::toDegrees;
        }
    }
}
