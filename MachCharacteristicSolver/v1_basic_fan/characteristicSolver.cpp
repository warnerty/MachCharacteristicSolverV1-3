#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

// Developed by **Tyler Warner**

constexpr double PI = 3.14159265358979323846;

// convert angles
double radians(double value) { return value * PI / 180.0; }
double degrees(double value) { return value * 180.0 / PI; }

// calculate the prandtl-meyer angle
double prandtlMeyer(double Mach, double gamma)
{
    if (!std::isfinite(Mach) || Mach < 1.0 ||
        !std::isfinite(gamma) || gamma <= 1.0)
        throw std::invalid_argument("Require Mach >= 1 and gamma > 1.");
    const double a = (gamma + 1.0) / (gamma - 1.0);
    const double b = Mach * Mach - 1.0;
    return std::sqrt(a) * std::atan(std::sqrt(b / a))
         - std::atan(std::sqrt(b));
}

// invert the angle using bisection
double machFromNu(double nu, double gamma)
{
    if (!std::isfinite(gamma) || gamma <= 1.0)
        throw std::invalid_argument("Require gamma > 1.");
    const double limit = PI / 2.0 *
        (std::sqrt((gamma + 1.0) / (gamma - 1.0)) - 1.0);
    if (!std::isfinite(nu) || nu < 0.0 || nu >= limit)
        throw std::invalid_argument("Expansion angle is outside the finite-Mach range.");
    if (nu == 0.0) return 1.0;

    double lower = 1.0, upper = 2.0;
    while (prandtlMeyer(upper, gamma) < nu)
    {
        upper *= 2.0;
        if (!std::isfinite(upper * upper))
            throw std::runtime_error("Unable to bracket Mach.");
    }
    for (int iteration = 0; iteration < 100; ++iteration)
    {
        const double middle = lower + (upper - lower) / 2.0;
        if (middle == lower || middle == upper) break;
        if (prandtlMeyer(middle, gamma) < nu) lower = middle;
        else upper = middle;
    }
    return lower + (upper - lower) / 2.0;
}

struct FanState
{
    double Mach, theta, nu, mu, alpha;
};

int main()
{
    try
    {
        // flow inputs
        const double gamma = 1.4;
        const double inletMach = 2.0;
        const double turnDegrees = 25.0;

        // geometry inputs
        const double cornerX = 0.0;
        const double cornerY = 0.0;
        const double rayLength = 10.0;
        const int intervals = 20;
        const int segmentsPerRay = 50;

        // check inputs before calculating
        if (intervals <= 0 || segmentsPerRay <= 0 ||
            !std::isfinite(rayLength) || rayLength <= 0.0 ||
            !std::isfinite(turnDegrees) || turnDegrees < 0.0 ||
            !std::isfinite(cornerX) || !std::isfinite(cornerY))
            throw std::invalid_argument("Invalid fan geometry.");

        const double totalTurn = radians(turnDegrees);
        const double inletNu = prandtlMeyer(inletMach, gamma);
        const double outletMach = machFromNu(inletNu + totalTurn, gamma);

        // calculate one characteristic
        auto calculateState = [&](int i)
        {
            const double delta =
                totalTurn * (static_cast<double>(i) / intervals);
            FanState state{};
            state.theta = -delta;
            state.nu = inletNu + delta;
            state.Mach = machFromNu(state.nu, gamma);
            state.mu = std::asin(1.0 / state.Mach);
            state.alpha = state.theta + state.mu;
            return state;
        };

        // save in the current working directory
        const auto outputPath = std::filesystem::absolute("characteristics.csv");
        std::ofstream csv(outputPath);
        if (!csv) throw std::runtime_error("Cannot open CSV output.");
        csv << std::setprecision(15);
        csv << "line_id,point_id,x,y,Mach,theta_deg,nu_deg,mu_deg,alpha_deg\n";

        // write rows in order
        for (int i = 0; i <= intervals; ++i)
        {
            const FanState state = calculateState(i);
            const double dx = std::cos(state.alpha);
            const double dy = std::sin(state.alpha);
            for (int j = 0; j <= segmentsPerRay; ++j)
            {
                const double s =
                    rayLength * (static_cast<double>(j) / segmentsPerRay);
                csv << i << ',' << j << ','
                    << cornerX + s * dx << ',' << cornerY + s * dy << ','
                    << state.Mach << ',' << degrees(state.theta) << ','
                    << degrees(state.nu) << ',' << degrees(state.mu) << ','
                    << degrees(state.alpha) << '\n';
            }
        }
        csv.close();
        if (!csv) throw std::runtime_error("CSV writing failed.");

        // print a short summary
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "Outlet Mach: " << outletMach << '\n';
        std::cout << "Saved: " << outputPath.string() << '\n';
    }
    catch (const std::exception& error)
    {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
