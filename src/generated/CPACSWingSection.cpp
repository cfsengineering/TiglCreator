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

#include <cassert>
#include "CCPACSWingSections.h"
#include "CPACSWingSection.h"
#include "CTiglError.h"
#include "CTiglLogging.h"
#include "CTiglUIDManager.h"
#include "TixiHelper.h"

namespace tigl
{
namespace generated
{
    CPACSWingSection::CPACSWingSection(CCPACSWingSections* parent, CTiglUIDManager* uidMgr)
        : m_uidMgr(uidMgr)
        , m_transformation(m_uidMgr)
        , m_elements(reinterpret_cast<CCPACSWingSection*>(this), m_uidMgr)
    {
        //assert(parent != NULL);
        m_parent = parent;
    }

    CPACSWingSection::~CPACSWingSection()
    {
        if (m_uidMgr) m_uidMgr->TryUnregisterObject(m_uID);
    }

    const CCPACSWingSections* CPACSWingSection::GetParent() const
    {
        return m_parent;
    }

    CCPACSWingSections* CPACSWingSection::GetParent()
    {
        return m_parent;
    }

    CTiglUIDManager& CPACSWingSection::GetUIDManager()
    {
        return *m_uidMgr;
    }

    const CTiglUIDManager& CPACSWingSection::GetUIDManager() const
    {
        return *m_uidMgr;
    }

    void CPACSWingSection::ReadCPACS(const TixiDocumentHandle& tixiHandle, const std::string& xpath)
    {
        // read attribute uID
        if (tixi::TixiCheckAttribute(tixiHandle, xpath, "uID")) {
            m_uID = tixi::TixiGetAttribute<std::string>(tixiHandle, xpath, "uID");
            if (m_uID.empty()) {
                LOG(WARNING) << "Required attribute uID is empty at xpath " << xpath;
            }
        }
        else {
            LOG(ERROR) << "Required attribute uID is missing at xpath " << xpath;
        }

        // read element name
        if (tixi::TixiCheckElement(tixiHandle, xpath + "/name")) {
            m_name = tixi::TixiGetElement<std::string>(tixiHandle, xpath + "/name");
            if (m_name.empty()) {
                LOG(WARNING) << "Required element name is empty at xpath " << xpath;
            }
        }
        else {
            LOG(ERROR) << "Required element name is missing at xpath " << xpath;
        }

        // read element description
        if (tixi::TixiCheckElement(tixiHandle, xpath + "/description")) {
            m_description = tixi::TixiGetElement<std::string>(tixiHandle, xpath + "/description");
            if (m_description->empty()) {
                LOG(WARNING) << "Optional element description is present but empty at xpath " << xpath;
            }
        }

        // read element transformation
        if (tixi::TixiCheckElement(tixiHandle, xpath + "/transformation")) {
            m_transformation.ReadCPACS(tixiHandle, xpath + "/transformation");
        }
        else {
            LOG(ERROR) << "Required element transformation is missing at xpath " << xpath;
        }

        // read element elements
        if (tixi::TixiCheckElement(tixiHandle, xpath + "/elements")) {
            m_elements.ReadCPACS(tixiHandle, xpath + "/elements");
        }
        else {
            LOG(ERROR) << "Required element elements is missing at xpath " << xpath;
        }

        if (m_uidMgr && !m_uID.empty()) m_uidMgr->RegisterObject(m_uID, *this);
    }

    void CPACSWingSection::WriteCPACS(const TixiDocumentHandle& tixiHandle, const std::string& xpath) const
    {
        // write attribute uID
        tixi::TixiSaveAttribute(tixiHandle, xpath, "uID", m_uID);

        // write element name
        tixi::TixiCreateElementIfNotExists(tixiHandle, xpath + "/name");
        tixi::TixiSaveElement(tixiHandle, xpath + "/name", m_name);

        // write element description
        if (m_description) {
            tixi::TixiCreateElementIfNotExists(tixiHandle, xpath + "/description");
            tixi::TixiSaveElement(tixiHandle, xpath + "/description", *m_description);
        }
        else {
            if (tixi::TixiCheckElement(tixiHandle, xpath + "/description")) {
                tixi::TixiRemoveElement(tixiHandle, xpath + "/description");
            }
        }

        // write element transformation
        tixi::TixiCreateElementIfNotExists(tixiHandle, xpath + "/transformation");
        m_transformation.WriteCPACS(tixiHandle, xpath + "/transformation");

        // write element elements
        tixi::TixiCreateElementIfNotExists(tixiHandle, xpath + "/elements");
        m_elements.WriteCPACS(tixiHandle, xpath + "/elements");

    }

    const std::string& CPACSWingSection::GetUID() const
    {
        return m_uID;
    }

    void CPACSWingSection::SetUID(const std::string& value)
    {
        if (m_uidMgr) {
            m_uidMgr->TryUnregisterObject(m_uID);
            m_uidMgr->RegisterObject(value, *this);
        }
        m_uID = value;
    }

    const std::string& CPACSWingSection::GetName() const
    {
        return m_name;
    }

    void CPACSWingSection::SetName(const std::string& value)
    {
        m_name = value;
    }

    const boost::optional<std::string>& CPACSWingSection::GetDescription() const
    {
        return m_description;
    }

    void CPACSWingSection::SetDescription(const boost::optional<std::string>& value)
    {
        m_description = value;
    }

    const CCPACSTransformation& CPACSWingSection::GetTransformation() const
    {
        return m_transformation;
    }

    CCPACSTransformation& CPACSWingSection::GetTransformation()
    {
        return m_transformation;
    }

    const CCPACSWingSectionElements& CPACSWingSection::GetElements() const
    {
        return m_elements;
    }

    CCPACSWingSectionElements& CPACSWingSection::GetElements()
    {
        return m_elements;
    }

} // namespace generated
} // namespace tigl
