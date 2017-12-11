// Copyright (c) 2016 RISC Software GmbH
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

#include "CPACSSparCell.h"
#include "CPACSSparCells.h"
#include "CTiglError.h"
#include "CTiglLogging.h"
#include "CTiglUIDManager.h"
#include "TixiHelper.h"

namespace tigl
{
    namespace generated
    {
        CPACSSparCells::CPACSSparCells(CTiglUIDManager* uidMgr) :
            m_uidMgr(uidMgr) {}
        
        CPACSSparCells::~CPACSSparCells() {}
        
        CTiglUIDManager& CPACSSparCells::GetUIDManager()
        {
            return *m_uidMgr;
        }
        
        const CTiglUIDManager& CPACSSparCells::GetUIDManager() const
        {
            return *m_uidMgr;
        }
        
        void CPACSSparCells::ReadCPACS(const TixiDocumentHandle& tixiHandle, const std::string& xpath)
        {
            // read element sparCell
            if (tixi::TixiCheckElement(tixiHandle, xpath + "/sparCell")) {
                tixi::TixiReadElements(tixiHandle, xpath + "/sparCell", m_sparCells, m_uidMgr);
            }
            
        }
        
        void CPACSSparCells::WriteCPACS(const TixiDocumentHandle& tixiHandle, const std::string& xpath) const
        {
            // write element sparCell
            tixi::TixiSaveElements(tixiHandle, xpath + "/sparCell", m_sparCells);
            
        }
        
        const std::vector<unique_ptr<CPACSSparCell> >& CPACSSparCells::GetSparCells() const
        {
            return m_sparCells;
        }
        
        std::vector<unique_ptr<CPACSSparCell> >& CPACSSparCells::GetSparCells()
        {
            return m_sparCells;
        }
        
        CPACSSparCell& CPACSSparCells::AddSparCell()
        {
            m_sparCells.push_back(make_unique<CPACSSparCell>(m_uidMgr));
            return *m_sparCells.back();
        }
        
        void CPACSSparCells::RemoveSparCell(CPACSSparCell& ref)
        {
            for (std::size_t i = 0; i < m_sparCells.size(); i++) {
                if (m_sparCells[i].get() == &ref) {
                    m_sparCells.erase(m_sparCells.begin() + i);
                    return;
                }
            }
            throw CTiglError("Element not found");
        }
        
    }
}
