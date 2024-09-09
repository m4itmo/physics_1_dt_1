#include <iostream>
#include <string>
#include <cmath>

using namespace std;

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
tuple<double, double, double> CartesianToSpherical(double x, double y, double z) {
    double r = sqrt(x*x + y*y + z*z);
    double theta = atan2(y, x);
    double phi = acos(z / r);

    return make_tuple(r, theta, phi);
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

string cout_tuple(tuple<double, double, double> t, string v1, string v2, string v3) {
    return v1 + " = " + to_string(get<0>(t)) + "\n" + v2 + " = " + to_string(get<1>(t)) + "\n" + v3 + " = " + to_string(get<2>(t));
}



int main() {
    cout << " ________  ________  ________  ________  ________  ___  ________   ________  _________" << endl;
    cout << "|\\   ____\\|\\   __  \\|\\   __  \\|\\   __  \\|\\   ___ \\|\\  \\|\\   ___  \\|\\   __  \\|\\___   ___\\" << endl;
    cout << "\\ \\  \\___|\\ \\  \\|\\  \\ \\  \\|\\  \\ \\  \\|\\  \\ \\  \\_|\\ \\ \\  \\ \\  \\\\ \\  \\ \\  \\|\\  \\|___ \\  \\_|" << endl;
    cout << " \\ \\  \\    \\ \\  \\\\\\  \\ \\  \\\\\\  \\ \\   _  _\\ \\  \\ \\\\ \\ \\  \\ \\  \\\\ \\  \\ \\   __  \\   \\ \\  \\" << endl;
    cout << "  \\ \\  \\____\\ \\  \\\\\\  \\ \\  \\\\\\  \\ \\  \\\\  \\\\ \\  \\_\\\\ \\ \\  \\ \\  \\\\ \\  \\ \\  \\ \\  \\   \\ \\  \\" << endl;
    cout << "   \\ \\_______\\ \\_______\\ \\_______\\ \\__\\\\ _\\\\ \\_______\\ \\__\\ \\__\\\\ \\__\\ \\__\\ \\__\\   \\ \\__\\" << endl;
    cout << "    \\|_______|\\|_______|\\|_______|\\|__|\\|__|\\|_______|\\|__|\\|__| \\|__|\\|__|\\|__|    \\|__|" << endl;
    cout << "" << endl << endl;
    cout << " _______           ________       ___    ___ ________  _________  _______   _____ ______" << endl;
    cout << "|\\  ___ \\         |\\   ____\\     |\\  \\  /  /|\\   ____\\|\\___   ___\\\\  ___ \\ |\\   _ \\  _   \\" << endl;
    cout << "\\ \\   __/|        \\ \\  \\___|_    \\ \\  \\/  / | \\  \\___|\\|___ \\  \\_\\ \\   __/|\\ \\  \\\\\\__\\ \\  \\" << endl;
    cout << " \\ \\  \\_|/__       \\ \\_____  \\    \\ \\    / / \\ \\_____  \\   \\ \\  \\ \\ \\  \\_|/_\\ \\  \\\\|__| \\  \\" << endl;
    cout << "  \\ \\  \\_|\\ \\       \\|____|\\  \\    \\/  /  /   \\|____|\\  \\   \\ \\  \\ \\ \\  \\_|\\ \\ \\  \\    \\ \\  \\" << endl;
    cout << "   \\ \\_______\\        ____\\_\\  \\ __/  / /       ____\\_\\  \\   \\ \\__\\ \\ \\_______\\ \\__\\    \\ \\__\\" << endl;
    cout << "    \\|_______|       |\\_________\\\\___/ /       |\\_________\\   \\|__|  \\|_______|\\|__|     \\|__|" << endl;
    cout << "                     \\|_________\\|___|/        \\|_________|" << endl;
    cout << "" << endl << endl;
    cout << " ________  ________  ________   ___      ___ _______   ________  _________  _______   ________" << endl;
    cout << "|\\   ____\\|\\   __  \\|\\   ___  \\|\\  \\    /  /|\\  ___ \\ |\\   __  \\|\\___   ___\\\\  ___ \\ |\\   __  \\" << endl;
    cout << "\\ \\  \\___|\\ \\  \\|\\  \\ \\  \\\\ \\  \\ \\  \\  /  / | \\   __/|\\ \\  \\|\\  \\|___ \\  \\_\\ \\   __/|\\ \\  \\|\\  \\" << endl;
    cout << " \\ \\  \\    \\ \\  \\\\\\  \\ \\  \\\\ \\  \\ \\  \\/  / / \\ \\  \\_|/_\\ \\   _  _\\   \\ \\  \\ \\ \\  \\_|/_\\ \\   _  _\\" << endl;
    cout << "  \\ \\  \\____\\ \\  \\\\\\  \\ \\  \\\\ \\  \\ \\    / /   \\ \\  \\_|\\ \\ \\  \\\\  \\|   \\ \\  \\ \\ \\  \\_|\\ \\ \\  \\\\  \\|" << endl;
    cout << "   \\ \\_______\\ \\_______\\ \\__\\\\ \\__\\ \\__/ /     \\ \\_______\\ \\__\\\\ _\\    \\ \\__\\ \\ \\_______\\ \\__\\\\ _\\" << endl;
    cout << "    \\|_______|\\|_______|\\|__| \\|__|\\|__|/       \\|_______|\\|__|\\|__|    \\|__|  \\|_______|\\|__|\\|__" << endl;
    cout << "" << endl << endl;
    cout << "github  :  https://github.com/m4itmo/physics_1_dt_1" << endl;
    cout << "author  :  std46d6b <dev@m5k.ru>" << endl;
    cout << "version :  0.1.1" << endl;
    cout << "" << endl << endl;


    int a, b, c;

    cout << "1 - Cartesian coordinate system -> Spherical coordinate system" << endl;
    cout << "2 - Spherical coordinate system -> Сartesian coordinate system" << endl;
    cout << "select option: ";
    
    int option;
    cin >> option;
    cout << endl;
    
    if (option == 1) {
        cout << "Cartesian coordinate system -> Spherical coordinate system" << endl;
        cout << "r = \\sqrt{x^2 + y^2 + z^2}" << endl;
        cout << "theta = arctg(\\sqrt{x^2 + y^2} / z)" << endl;
        cout << "phi = arctg(y / x)" << endl;
        cout << "" << endl;
        cout << "Enter x, y, z: ";
        cin >> a >> b >> c;
        
        tuple<double, double, double> result = CartesianToSpherical(a, b, c);
        cout << cout_tuple(result, "r    ", "theta", "phi  ") << endl;
    } else if (option == 2) {
        cout << "Spherical coordinate system -> Cartesian coordinate system" << endl;
        cout << "x = r * sin(theta) * cos(phi)" << endl;
        cout << "y = r * sin(theta) * sin(phi)" << endl;
        cout << "z = r * cos(theta)" << endl;
        cout << "" << endl;
        cout << "Enter r, theta, phi: ";
        cin >> a >> b >> c;

        tuple<double, double, double> result = SphericalToCartesian(a, b, c);
        cout << cout_tuple(result, "x    ", "y    ", "z    ") << endl;
    }

    return 0;
}
