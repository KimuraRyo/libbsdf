// =================================================================== //
// Copyright (C) 2019-2023 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Common/Array.h>

#include <libbsdf/Common/Utility.h>

using namespace lb;

void array_util::findBounds(const Arrayd& values,
                            double        value,
                            bool          equalIntervalValues,
                            int*          lowerIndex,
                            int*          upperIndex,
                            double*       lowerValue,
                            double*       upperValue)
{
    if (values.size() == 1) {
        *lowerIndex = 0;
        *upperIndex = 0;
        *lowerValue = values[0];
        *upperValue = values[0];

        return;
    }

    int backIndex = static_cast<int>(values.size() - 1);
    if (equalIntervalValues) {
        // Calculate lower and upper indices.
        *lowerIndex = static_cast<int>(backIndex * (value / values[backIndex]));
        *lowerIndex = std::min(*lowerIndex, backIndex - 1);
        *upperIndex = *lowerIndex + 1;
    }
    else {
        // Find lower and upper indices.
        const double* valuePtr = std::lower_bound(&values[0], &values[0] + values.size(), value);
        *upperIndex = clamp(static_cast<int>(valuePtr - &values[0]), 1, backIndex);
        *lowerIndex = *upperIndex - 1;
    }

    *lowerValue = values[*lowerIndex];
    *upperValue = values[*upperIndex];
}
