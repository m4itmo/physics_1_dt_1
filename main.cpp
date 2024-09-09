#include <iostream>
#include <string>
#include <cmath>
#include <tuple>

/**
 * Converts Cartesian coordinates to Spherical coordinates.
 *
 * @param x x-coordinate in Cartesian system
 * @param y y-coordinate in Cartesian system
 * @param z z-coordinate in Cartesian system
 * @return tuple containing (r, theta, phi) in Spherical coordinate system
 *
 * Formulas:
 *   r = sqrt(x**2 + y**2 + z**2)
 *   theta = atan2(y, x)
 *   phi = acos(z / r)
 */
std::tuple<double, double, double> CartesianToSpherical(double x, double y, double z) {
    double r = std::sqrt(x * x + y * y + z * z);
    double theta = std::atan2(y, x);
    double phi = std::acos(z / r);

    return std::make_tuple(r, theta, phi);
}

/**
 * Converts Spherical coordinates to Cartesian coordinates.
 *
 * @param r radial distance from the origin
 * @param theta azimuthal angle in the xy-plane from the positive x-axis
 * @param phi polar angle from the positive z-axis
 * @return std::tuple containing (x, y, z) in Cartesian coordinate system
 *
 * Formulas:
 *   x = r * sin(phi) * cos(theta)
 *   y = r * sin(phi) * sin(theta)
 *   z = r * cos(phi)
 */
std::tuple<double, double, double> SphericalToCartesian(double r, double theta, double phi) {
    double x = r * std::sin(phi) * std::cos(theta);
    double y = r * std::sin(phi) * std::sin(theta);
    double z = r * std::cos(phi);

    return std::make_tuple(x, y, z);
}

/**
 * Converts Cylindrical coordinates to Cartesian coordinates.
 *
 * @param r radial distance from the origin (positive)
 * @param theta azimuthal angle in the xy-plane from the positive x-axis
 * @param z z-coordinate in Cartesian system
 * @return std::tuple containing (x, y, z) in Cartesian coordinate system
 *
 * Formulas:
 *   x = r * cos(theta)
 *   y = r * sin(theta)
 *   z = z
 */
std::tuple<double, double, double> CylindricalToCartesian(double r, double theta, double z) {
    double x = r * std::cos(theta);
    double y = r * std::sin(theta);

    return std::make_tuple(x, y, z);
}

/**
 * Converts Cartesian coordinates to Cylindrical coordinates.
 *
 * @param x x-coordinate in Cartesian system
 * @param y y-coordinate in Cartesian system
 * @param z z-coordinate in Cartesian system
 * @return std::tuple containing (r, theta, z) in Cylindrical coordinate system
 *
 * Formulas:
 *   r = sqrt(x**2 + y**2)
 *   theta = atan2(y, x)
 *   z = z
 */
std::tuple<double, double, double> CartesianToCylindrical(double x, double y, double z) {
    double r = std::sqrt(x * x + y * y);
    double theta = std::atan2(y, x);

    return std::make_tuple(r, theta, z);
}

/**
 * Converts Cylindrical coordinates to Spherical coordinates.
 *
 * @param r radial distance from the origin (positive)
 * @param theta azimuthal angle in the xy-plane from the positive x-axis
 * @param z z-coordinate in Cartesian system
 * @return std::tuple containing (r, theta, phi) in Spherical coordinate system
 *
 * Formulas:
 *   r = sqrt(x**2 + y**2 + z**2)
 *   theta = atan2(y, x)
 *   phi = acos(z / r)
 */
std::tuple<double, double, double> CylindricalToSpherical(double r, double theta, double z) {
    double radius = std::sqrt(r * r + z * z);
    double phi = std::acos(z / radius);

    return std::make_tuple(radius, theta, phi);
}

/**
 * Converts Spherical coordinates to Cylindrical coordinates.
 *
 * @param r radial distance from the origin (positive)
 * @param theta azimuthal angle in the xy-plane from the positive x-axis
 * @param phi polar angle from the positive z-axis
 * @return std::tuple containing (r, theta, z) in Cylindrical coordinate system
 *
 * Formulas:
 *   r = r * sin(phi)
 *   theta = theta
 *   z = r * cos(phi)
 */
std::tuple<double, double, double> SphericalToCylindrical(double r, double theta, double phi) {
    double radial = r * std::sin(phi);

    return std::make_tuple(radial, theta, r * std::cos(phi));
}

/**
 * Converts tuple to string.
 *
 * @param t tuple
 * @param v1 first variable
 * @param v2 second variable
 * @param v3 third variable
 * @return string
 */
std::string cout_tuple(std::tuple<double, double, double> t, std::string v1, std::string v2, std::string v3) {
    return v1 + " = " + std::to_string(std::get<0>(t)) + "\n" + v2 + " = " + std::to_string(std::get<1>(t)) + "\n" + v3 + " = " +
           std::to_string(std::get<2>(t));
}

int main() {
    std::cout << " ________  ________  ________  ________  ________  ___  ________   ________  _________" << std::endl;
    std::cout << "|\\   ____\\|\\   __  \\|\\   __  \\|\\   __  \\|\\   ___ \\|\\  \\|\\   ___  \\|\\   __  \\|\\___   ___\\" << std::endl;
    std::cout << "\\ \\  \\___|\\ \\  \\|\\  \\ \\  \\|\\  \\ \\  \\|\\  \\ \\  \\_|\\ \\ \\  \\ \\  \\\\ \\  \\ \\  \\|\\  \\|___ \\  \\_|" << std::endl;
    std::cout << " \\ \\  \\    \\ \\  \\\\\\  \\ \\  \\\\\\  \\ \\   _  _\\ \\  \\ \\\\ \\ \\  \\ \\  \\\\ \\  \\ \\   __  \\   \\ \\  \\" << std::endl;
    std::cout << "  \\ \\  \\____\\ \\  \\\\\\  \\ \\  \\\\\\  \\ \\  \\\\  \\\\ \\  \\_\\\\ \\ \\  \\ \\  \\\\ \\  \\ \\  \\ \\  \\   \\ \\  \\" << std::endl;
    std::cout << "   \\ \\_______\\ \\_______\\ \\_______\\ \\__\\\\ _\\\\ \\_______\\ \\__\\ \\__\\\\ \\__\\ \\__\\ \\__\\   \\ \\__\\" << std::endl;
    std::cout << "    \\|_______|\\|_______|\\|_______|\\|__|\\|__|\\|_______|\\|__|\\|__| \\|__|\\|__|\\|__|    \\|__|" << std::endl;
    std::cout << "" << std::endl << std::endl;
    std::cout << " _______           ________       ___    ___ ________  _________  _______   _____ ______" << std::endl;
    std::cout << "|\\  ___ \\         |\\   ____\\     |\\  \\  /  /|\\   ____\\|\\___   ___\\\\  ___ \\ |\\   _ \\  _   \\" << std::endl;
    std::cout << "\\ \\   __/|        \\ \\  \\___|_    \\ \\  \\/  / | \\  \\___|\\|___ \\  \\_\\ \\   __/|\\ \\  \\\\\\__\\ \\  \\" << std::endl;
    std::cout << " \\ \\  \\_|/__       \\ \\_____  \\    \\ \\    / / \\ \\_____  \\   \\ \\  \\ \\ \\  \\_|/_\\ \\  \\\\|__| \\  \\" << std::endl;
    std::cout << "  \\ \\  \\_|\\ \\       \\|____|\\  \\    \\/  /  /   \\|____|\\  \\   \\ \\  \\ \\ \\  \\_|\\ \\ \\  \\    \\ \\  \\" << std::endl;
    std::cout << "   \\ \\_______\\        ____\\_\\  \\ __/  / /       ____\\_\\  \\   \\ \\__\\ \\ \\_______\\ \\__\\    \\ \\__\\" << std::endl;
    std::cout << "    \\|_______|       |\\_________\\\\___/ /       |\\_________\\   \\|__|  \\|_______|\\|__|     \\|__|" << std::endl;
    std::cout << "                     \\|_________\\|___|/        \\|_________|" << std::endl;
    std::cout << "" << std::endl << std::endl;
    std::cout << " ________  ________  ________   ___      ___ _______   ________  _________  _______   ________" << std::endl;
    std::cout << "|\\   ____\\|\\   __  \\|\\   ___  \\|\\  \\    /  /|\\  ___ \\ |\\   __  \\|\\___   ___\\\\  ___ \\ |\\   __  \\" << std::endl;
    std::cout << "\\ \\  \\___|\\ \\  \\|\\  \\ \\  \\\\ \\  \\ \\  \\  /  / | \\   __/|\\ \\  \\|\\  \\|___ \\  \\_\\ \\   __/|\\ \\  \\|\\  \\" << std::endl;
    std::cout << " \\ \\  \\    \\ \\  \\\\\\  \\ \\  \\\\ \\  \\ \\  \\/  / / \\ \\  \\_|/_\\ \\   _  _\\   \\ \\  \\ \\ \\  \\_|/_\\ \\   _  _\\" << std::endl;
    std::cout << "  \\ \\  \\____\\ \\  \\\\\\  \\ \\  \\\\ \\  \\ \\    / /   \\ \\  \\_|\\ \\ \\  \\\\  \\|   \\ \\  \\ \\ \\  \\_|\\ \\ \\  \\\\  \\|" << std::endl;
    std::cout << "   \\ \\_______\\ \\_______\\ \\__\\\\ \\__\\ \\__/ /     \\ \\_______\\ \\__\\\\ _\\    \\ \\__\\ \\ \\_______\\ \\__\\\\ _\\" << std::endl;
    std::cout << "    \\|_______|\\|_______|\\|__| \\|__|\\|__|/       \\|_______|\\|__|\\|__|    \\|__|  \\|_______|\\|__|\\|__" << std::endl;
    std::cout << "" << std::endl << std::endl;
    std::cout << "github  :  https://github.com/m4itmo/physics_1_dt_1" << std::endl;
    std::cout << "author  :  std46d6b <dev@m5k.ru>" << std::endl;
    std::cout << "version :  0.2.0" << std::endl;
    std::cout << "" << std::endl << std::endl;


    int a, b, c;

    std::cout << "1 - Cartesian coordinate system   -> Spherical coordinate system" << std::endl;
    std::cout << "2 - Spherical coordinate system   -> Сartesian coordinate system" << std::endl;

    std::cout << "3 - Cylindrical coordinate system -> Cartesian coordinate system" << std::endl;
    std::cout << "4 - Cartesian coordinate system   -> Cylindrical coordinate system" << std::endl;
    std::cout << "5 - Cylindrical coordinate system -> Spherical coordinate system" << std::endl;
    std::cout << "6 - Spherical coordinate system   -> Cylindrical coordinate system" << std::endl;
    std::cout << "select option: ";

    int option;
    std::cin >> option;
    std::cout << std::endl;

    if (option == 1) {
        std::cout << "Cartesian coordinate system -> Spherical coordinate system" << std::endl;
        std::cout << "r = \\sqrt{x^2 + y^2 + z^2}" << std::endl;
        std::cout << "theta = arctg(\\sqrt{x^2 + y^2} / z)" << std::endl;
        std::cout << "phi = arctg(y / x)" << std::endl;
        std::cout << "" << std::endl;
        std::cout << "Enter x, y, z: ";
        std::cin >> a >> b >> c;

        std::tuple<double, double, double> result = CartesianToSpherical(a, b, c);
        std::cout << cout_tuple(result, "r    ", "theta", "phi  ") << std::endl;
    } else if (option == 2) {
        std::cout << "Spherical coordinate system -> Cartesian coordinate system" << std::endl;
        std::cout << "x = r * sin(theta) * cos(phi)" << std::endl;
        std::cout << "y = r * sin(theta) * sin(phi)" << std::endl;
        std::cout << "z = r * cos(theta)" << std::endl;
        std::cout << "" << std::endl;
        std::cout << "Enter r, theta, phi: ";
        std::cin >> a >> b >> c;

        std::tuple<double, double, double> result = SphericalToCartesian(a, b, c);
        std::cout << cout_tuple(result, "x    ", "y    ", "z    ") << std::endl;
    } else if (option == 3) {
        std::cout << "Cylindrical coordinate system -> Cartesian coordinate system" << std::endl;
        std::cout << "x = r * cos(phi)" << std::endl;
        std::cout << "y = r * sin(phi)" << std::endl;
        std::cout << "z = z" << std::endl;
        std::cout << "" << std::endl;
        std::cout << "Enter r, phi, z: ";
        std::cin >> a >> b >> c;

        std::tuple<double, double, double> result = CylindricalToCartesian(a, b, c);
        std::cout << cout_tuple(result, "x    ", "y    ", "z    ") << std::endl;
    } else if (option == 4) {
        std::cout << "Cartesian coordinate system -> Cylindrical coordinate system" << std::endl;
        std::cout << "r = \\sqrt{x^2 + y^2}" << std::endl;
        std::cout << "phi = arctg(y / x)" << std::endl;
        std::cout << "z = z" << std::endl;
        std::cout << "" << std::endl;
        std::cout << "Enter x, y, z: ";
        std::cin >> a >> b >> c;

        std::tuple<double, double, double> result = CartesianToCylindrical(a, b, c);
        std::cout << cout_tuple(result, "r    ", "phi  ", "z    ") << std::endl;
    } else if (option == 5) {
        std::cout << "Cylindrical coordinate system -> Spherical coordinate system" << std::endl;
        std::cout << "r = \\sqrt{x^2 + y^2}" << std::endl;
        std::cout << "theta = arctg(\\sqrt{x^2 + y^2} / z)" << std::endl;
        std::cout << "phi = arctg(y / x)" << std::endl;
        std::cout << "" << std::endl;
        std::cout << "Enter r, phi, z: ";
        std::cin >> a >> b >> c;

        std::tuple<double, double, double> result = CylindricalToSpherical(a, b, c);
        std::cout << cout_tuple(result, "r    ", "theta", "phi  ") << std::endl;
    } else if (option == 6) {
        std::cout << "Spherical coordinate system -> Cylindrical coordinate system" << std::endl;
        std::cout << "x = r * sin(theta) * cos(phi)" << std::endl;
        std::cout << "y = r * sin(theta) * sin(phi)" << std::endl;
        std::cout << "z = r * cos(theta)" << std::endl;
        std::cout << "" << std::endl;
        std::cout << "Enter r, theta, phi: ";
        std::cin >> a >> b >> c;

        std::tuple<double, double, double> result = SphericalToCylindrical(a, b, c);
        std::cout << cout_tuple(result, "x    ", "y    ", "z    ") << std::endl;
    }
    system("pause");
    return 0;
}