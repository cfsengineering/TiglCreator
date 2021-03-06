// Copyright (c) 2018 RISC Software GmbH
//
// This file was generated by CPACSGen from CPACS XML Schema (c) German Aerospace Center (DLR/SC).
// Do not edit, all changes are lost when files are re-generated.
//
// Licensed under the Apache License, Version 2.0 (the "License")
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "CPACSWingRibPoint.h"
#include "CTiglError.h"
#include "CTiglLogging.h"
#include "TixiHelper.h"

namespace tigl
{
namespace generated
{
    CPACSWingRibPoint::CPACSWingRibPoint()
        : m_xsi(0)
    {
    }

    CPACSWingRibPoint::~CPACSWingRibPoint()
    {
    }

    void CPACSWingRibPoint::ReadCPACS(const TixiDocumentHandle& tixiHandle, const std::string& xpath)
    {
        // read element ribDefinitionUID
        if (tixi::TixiCheckElement(tixiHandle, xpath + "/ribDefinitionUID")) {
            m_ribDefinitionUID = tixi::TixiGetElement<std::string>(tixiHandle, xpath + "/ribDefinitionUID");
            if (m_ribDefinitionUID.empty()) {
                LOG(WARNING) << "Required element ribDefinitionUID is empty at xpath " << xpath;
            }
        }
        else {
            LOG(ERROR) << "Required element ribDefinitionUID is missing at xpath " << xpath;
        }

        // read element ribNumber
        if (tixi::TixiCheckElement(tixiHandle, xpath + "/ribNumber")) {
            m_ribNumber = tixi::TixiGetElement<int>(tixiHandle, xpath + "/ribNumber");
        }

        // read element xsi
        if (tixi::TixiCheckElement(tixiHandle, xpath + "/xsi")) {
            m_xsi = tixi::TixiGetElement<double>(tixiHandle, xpath + "/xsi");
        }
        else {
            LOG(ERROR) << "Required element xsi is missing at xpath " << xpath;
        }

    }

    void CPACSWingRibPoint::WriteCPACS(const TixiDocumentHandle& tixiHandle, const std::string& xpath) const
    {
        // write element ribDefinitionUID
        tixi::TixiCreateElementIfNotExists(tixiHandle, xpath + "/ribDefinitionUID");
        tixi::TixiSaveElement(tixiHandle, xpath + "/ribDefinitionUID", m_ribDefinitionUID);

        // write element ribNumber
        if (m_ribNumber) {
            tixi::TixiCreateElementIfNotExists(tixiHandle, xpath + "/ribNumber");
            tixi::TixiSaveElement(tixiHandle, xpath + "/ribNumber", *m_ribNumber);
        }
        else {
            if (tixi::TixiCheckElement(tixiHandle, xpath + "/ribNumber")) {
                tixi::TixiRemoveElement(tixiHandle, xpath + "/ribNumber");
            }
        }

        // write element xsi
        tixi::TixiCreateElementIfNotExists(tixiHandle, xpath + "/xsi");
        tixi::TixiSaveElement(tixiHandle, xpath + "/xsi", m_xsi);

    }

    const std::string& CPACSWingRibPoint::GetRibDefinitionUID() const
    {
        return m_ribDefinitionUID;
    }

    void CPACSWingRibPoint::SetRibDefinitionUID(const std::string& value)
    {
        m_ribDefinitionUID = value;
    }

    const boost::optional<int>& CPACSWingRibPoint::GetRibNumber() const
    {
        return m_ribNumber;
    }

    void CPACSWingRibPoint::SetRibNumber(const boost::optional<int>& value)
    {
        m_ribNumber = value;
    }

    const double& CPACSWingRibPoint::GetXsi() const
    {
        return m_xsi;
    }

    void CPACSWingRibPoint::SetXsi(const double& value)
    {
        m_xsi = value;
    }

} // namespace generated
} // namespace tigl
