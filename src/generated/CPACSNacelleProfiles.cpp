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

#include <CCPACSNacelleProfile.h>
#include "CPACSNacelleProfiles.h"
#include "CTiglError.h"
#include "CTiglLogging.h"
#include "CTiglUIDManager.h"
#include "TixiHelper.h"

namespace tigl
{
namespace generated
{
    CPACSNacelleProfiles::CPACSNacelleProfiles(CTiglUIDManager* uidMgr)
        : m_uidMgr(uidMgr)
    {
    }

    CPACSNacelleProfiles::~CPACSNacelleProfiles()
    {
    }

    CTiglUIDManager& CPACSNacelleProfiles::GetUIDManager()
    {
        return *m_uidMgr;
    }

    const CTiglUIDManager& CPACSNacelleProfiles::GetUIDManager() const
    {
        return *m_uidMgr;
    }

    void CPACSNacelleProfiles::ReadCPACS(const TixiDocumentHandle& tixiHandle, const std::string& xpath)
    {
        // read element nacelleProfile
        if (tixi::TixiCheckElement(tixiHandle, xpath + "/nacelleProfile")) {
            tixi::TixiReadElements(tixiHandle, xpath + "/nacelleProfile", m_nacelleProfiles, m_uidMgr);
        }

    }

    void CPACSNacelleProfiles::WriteCPACS(const TixiDocumentHandle& tixiHandle, const std::string& xpath) const
    {
        // write element nacelleProfile
        tixi::TixiSaveElements(tixiHandle, xpath + "/nacelleProfile", m_nacelleProfiles);

    }

    const std::vector<unique_ptr<CCPACSNacelleProfile> >& CPACSNacelleProfiles::GetNacelleProfiles() const
    {
        return m_nacelleProfiles;
    }

    std::vector<unique_ptr<CCPACSNacelleProfile> >& CPACSNacelleProfiles::GetNacelleProfiles()
    {
        return m_nacelleProfiles;
    }

    CCPACSNacelleProfile& CPACSNacelleProfiles::AddNacelleProfile()
    {
        m_nacelleProfiles.push_back(make_unique<CCPACSNacelleProfile>(m_uidMgr));
        return *m_nacelleProfiles.back();
    }

    void CPACSNacelleProfiles::RemoveNacelleProfile(CCPACSNacelleProfile& ref)
    {
        for (std::size_t i = 0; i < m_nacelleProfiles.size(); i++) {
            if (m_nacelleProfiles[i].get() == &ref) {
                m_nacelleProfiles.erase(m_nacelleProfiles.begin() + i);
                return;
            }
        }
        throw CTiglError("Element not found");
    }

} // namespace generated
} // namespace tigl
