/*

8888888b.   .d8888b.  8888888b.
888   Y88b d88P  Y88b 888   Y88b
888    888 888    888 888    888
888   d88P 888        888   d88P 888d888 .d88b.  88888b.  .d8888b
8888888P"  888        8888888P"  888P"  d88""88b 888 "88b 88K
888        888    888 888        888    888  888 888  888 "Y8888b.
888        Y88b  d88P 888        888    Y88..88P 888 d88P      X88
888         "Y8888P"  888        888     "Y88P"  88888P"   88888P'
                                                 888
                                                 888
                                                 888

Copyright (c) 2020 Kenneth Troldal Balslev

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

#ifndef PCPROPS_PENGROBINSON_HPP
#define PCPROPS_PENGROBINSON_HPP

#include <functional>
#include <memory>
#include <tuple>
#include <string>

namespace PCProps::EquationOfState
{
    /**
     * @brief
     */
    class PengRobinson
    {
        using JSONString = std::string;

    public:
        // =====================================================================
        // CONSTRUCTORS, ASSIGNMENT & INITIATION
        // =====================================================================

        /**
         * @brief Default constructor. The object only supports copying and moving. The init
         * function must be called before using the other member functions.
         */
        PengRobinson();

        /**
         * @brief
         * @param pureComponent
         */
        PengRobinson(const std::function<double(std::string)>& constants, const std::function<double(std::string, double)>& correlations);

        /**
         * @brief Copy constructor
         */
        PengRobinson(const PengRobinson& other);

        /**
         * @brief Move constructor
         */
        PengRobinson(PengRobinson&& other) noexcept;

        /**
         * @brief Destructor
         */
        ~PengRobinson();

        /**
         * @brief Copy assignment operator
         */
        PengRobinson& operator=(const PengRobinson& other);

        /**
         * @brief Move assignment operator
         */
        PengRobinson& operator=(PengRobinson&& other) noexcept;

        /**
         *
         * @tparam PC
         * @param pureComponent
         */
        template<typename PC>
        void init(const PC& pureComponent) {
            init([&](const std::string& ID)->double {return pureComponent.property(ID);},
                 [&](const std::string& ID, double t)->double {return pureComponent.correlation(ID, t);});
        }


        // =====================================================================
        // FLASH ALGORITHMS
        // =====================================================================

        /**
         * @brief
         * @param specification
         * @param var1
         * @param var2
         * @return
         */
        JSONString flash(const std::string& specification, double var1, double var2) const;

        /**
         * @brief Calculate the saturation pressure at the given temperature.
         * @param temperature The temperature [K]
         * @return The saturation pressure [Pa]
         * @warning Returns NaN if the temperature is higher than the critical temperature.
         */
        double computePSat(double temperature) const;

        /**
         * @brief Calculate the saturation temperature at the given pressure.
         * @param pressure The pressure [Pa]
         * @return The saturation temperature [K]
         * @warning Returns NaN if pressure is higher than the critical pressure.
         */
        double computeTSat(double pressure) const;

        JSONString computeProperties(double pressure, double temperature) const;


    private:

        /**
         * @brief Initiates an existing object with a new pure component.
         * @param pureComponent An object with an interface compatible with IPureComponent
         */
        void init(const std::function<double(std::string)>& constants, const std::function<double(std::string, double)>& correlations);

        class impl;
        std::unique_ptr<impl> m_impl;
    };

}    // namespace PCProps::PropertyPackage
#endif    // PCPROPS_PENGROBINSON_HPP
