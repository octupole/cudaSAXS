#include "Cell.h"
#include <cmath>
#include <stdexcept>
#include "Ftypedefs.h"

std::vector<std::vector<float>> Cell::co = std::vector<std::vector<float>>(DIM, std::vector<float>(DIM, 0.0f));
std::vector<std::vector<float>> Cell::oc = std::vector<std::vector<float>>(DIM, std::vector<float>(DIM, 0.0f));
void Cell::calculateMatrices(
    float a, float b, float c, float alpha, float beta, float gamma)
{
    // Error handling for invalid cell parameters
    if (a <= 0 || b <= 0 || c <= 0)
    {
        throw std::invalid_argument("Cell parameters must be positive.");
    }
    if (alpha <= 0 || alpha >= 180 || beta <= 0 || beta >= 180 || gamma <= 0 || gamma >= 180)
    {
        throw std::invalid_argument("Angles must be between 0 and 180 degrees.");
    }

    // Convert angles to radians
    float alphaRad = alpha * M_PI / 180.0;
    float betaRad = beta * M_PI / 180.0;
    float gammaRad = gamma * M_PI / 180.0;

    // Calculate cosines and sines
    float cosAlpha = std::cos(alphaRad);
    float cosBeta = std::cos(betaRad);
    float cosGamma = std::cos(gammaRad);
    float sinGamma = std::sin(gammaRad);

    // Transformation matrix (co)
    co = {
        {a, 0.0, 0.0},
        {b * cosGamma, b * sinGamma, 0.0},
        {c * cosBeta, c * (cosAlpha - cosBeta * cosGamma) / sinGamma, c * std::sqrt(1 - cosAlpha * cosAlpha - cosBeta * cosBeta - cosGamma * cosGamma + 2 * cosAlpha * cosBeta * cosGamma) / sinGamma}};

    // Calculate inverse matrix elements (oc)
    float V = co[0][0] * co[1][1] * co[2][2]; // Volume of the unit cell
    oc = {
        {1.0f / a, 0.0f, 0.0f},
        {(-cosGamma) / (a * sinGamma), 1.0f / (b * sinGamma), 0.0f},
        {b * c * (cosAlpha * cosGamma - cosBeta) / V, a * c * (cosBeta * cosGamma - cosAlpha) / V, a * b * sinGamma / V}};
}
