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

#include "CPACSMaterialDefinitionForProfileBased.h"
#include "CPACSMaterialDefinitionForProfileBasedPoint.h"
#include "CPACSProfileBasedStructuralElement.h"
#include "CTiglError.h"
#include "CTiglLogging.h"
#include "CTiglUIDManager.h"
#include "TixiHelper.h"

namespace tigl
{
namespace generated
{
    CPACSProfileBasedStructuralElement::CPACSProfileBasedStructuralElement(CTiglUIDManager* uidMgr)
        : m_uidMgr(uidMgr)
    {
    }

    CPACSProfileBasedStructuralElement::~CPACSProfileBasedStructuralElement()
    {
        if (m_uidMgr) m_uidMgr->TryUnregisterObject(m_uID);
    }

    CTiglUIDManager& CPACSProfileBasedStructuralElement::GetUIDManager()
    {
        return *m_uidMgr;
    }

    const CTiglUIDManager& CPACSProfileBasedStructuralElement::GetUIDManager() const
    {
        return *m_uidMgr;
    }

    void CPACSProfileBasedStructuralElement::ReadCPACS(const TixiDocumentHandle& tixiHandle, const std::string& xpath)
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
            if (m_name->empty()) {
                LOG(WARNING) << "Optional element name is present but empty at xpath " << xpath;
            }
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
            m_transformation = boost::in_place(m_uidMgr);
            try {
                m_transformation->ReadCPACS(tixiHandle, xpath + "/transformation");
            } catch(const std::exception& e) {
                LOG(ERROR) << "Failed to read transformation at xpath " << xpath << ": " << e.what();
                m_transformation = boost::none;
            }
        }

        // read element globalBeamProperties
        if (tixi::TixiCheckElement(tixiHandle, xpath + "/globalBeamProperties")) {
            m_globalBeamProperties_choice1 = boost::in_place(m_uidMgr);
            try {
                m_globalBeamProperties_choice1->ReadCPACS(tixiHandle, xpath + "/globalBeamProperties");
            } catch(const std::exception& e) {
                LOG(ERROR) << "Failed to read globalBeamProperties at xpath " << xpath << ": " << e.what();
                m_globalBeamProperties_choice1 = boost::none;
            }
        }

        // read element sheetProperties
        if (tixi::TixiCheckElement(tixiHandle, xpath + "/sheetProperties")) {
            tixi::TixiReadElements(tixiHandle, xpath + "/sheetProperties", m_sheetProperties_choice2);
        }

        // read element standardProfileType
        if (tixi::TixiCheckElement(tixiHandle, xpath + "/standardProfileType")) {
            m_standardProfileType_choice2_1 = stringToCPACSProfileBasedStructuralElement_standardProfileType(tixi::TixiGetElement<std::string>(tixiHandle, xpath + "/standardProfileType"));
        }

        // read element structuralProfileUID
        if (tixi::TixiCheckElement(tixiHandle, xpath + "/structuralProfileUID")) {
            m_structuralProfileUID_choice2_2 = tixi::TixiGetElement<std::string>(tixiHandle, xpath + "/structuralProfileUID");
            if (m_structuralProfileUID_choice2_2->empty()) {
                LOG(WARNING) << "Optional element structuralProfileUID is present but empty at xpath " << xpath;
            }
        }

        // read element pointProperties
        if (tixi::TixiCheckElement(tixiHandle, xpath + "/pointProperties")) {
            tixi::TixiReadElements(tixiHandle, xpath + "/pointProperties", m_pointProperties_choice2_2);
        }

        // read element referencePointUID
        if (tixi::TixiCheckElement(tixiHandle, xpath + "/referencePointUID")) {
            m_referencePointUID_choice2_2 = tixi::TixiGetElement<std::string>(tixiHandle, xpath + "/referencePointUID");
            if (m_referencePointUID_choice2_2->empty()) {
                LOG(WARNING) << "Optional element referencePointUID is present but empty at xpath " << xpath;
            }
        }

        if (m_uidMgr && !m_uID.empty()) m_uidMgr->RegisterObject(m_uID, *this);
        if (!ValidateChoices()) {
            LOG(ERROR) << "Invalid choice configuration at xpath " << xpath;
        }
    }

    void CPACSProfileBasedStructuralElement::WriteCPACS(const TixiDocumentHandle& tixiHandle, const std::string& xpath) const
    {
        // write attribute uID
        tixi::TixiSaveAttribute(tixiHandle, xpath, "uID", m_uID);

        // write element name
        if (m_name) {
            tixi::TixiCreateElementIfNotExists(tixiHandle, xpath + "/name");
            tixi::TixiSaveElement(tixiHandle, xpath + "/name", *m_name);
        }
        else {
            if (tixi::TixiCheckElement(tixiHandle, xpath + "/name")) {
                tixi::TixiRemoveElement(tixiHandle, xpath + "/name");
            }
        }

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
        if (m_transformation) {
            tixi::TixiCreateElementIfNotExists(tixiHandle, xpath + "/transformation");
            m_transformation->WriteCPACS(tixiHandle, xpath + "/transformation");
        }
        else {
            if (tixi::TixiCheckElement(tixiHandle, xpath + "/transformation")) {
                tixi::TixiRemoveElement(tixiHandle, xpath + "/transformation");
            }
        }

        // write element globalBeamProperties
        if (m_globalBeamProperties_choice1) {
            tixi::TixiCreateElementIfNotExists(tixiHandle, xpath + "/globalBeamProperties");
            m_globalBeamProperties_choice1->WriteCPACS(tixiHandle, xpath + "/globalBeamProperties");
        }
        else {
            if (tixi::TixiCheckElement(tixiHandle, xpath + "/globalBeamProperties")) {
                tixi::TixiRemoveElement(tixiHandle, xpath + "/globalBeamProperties");
            }
        }

        // write element sheetProperties
        tixi::TixiSaveElements(tixiHandle, xpath + "/sheetProperties", m_sheetProperties_choice2);

        // write element standardProfileType
        if (m_standardProfileType_choice2_1) {
            tixi::TixiCreateElementIfNotExists(tixiHandle, xpath + "/standardProfileType");
            tixi::TixiSaveElement(tixiHandle, xpath + "/standardProfileType", CPACSProfileBasedStructuralElement_standardProfileTypeToString(*m_standardProfileType_choice2_1));
        }
        else {
            if (tixi::TixiCheckElement(tixiHandle, xpath + "/standardProfileType")) {
                tixi::TixiRemoveElement(tixiHandle, xpath + "/standardProfileType");
            }
        }

        // write element structuralProfileUID
        if (m_structuralProfileUID_choice2_2) {
            tixi::TixiCreateElementIfNotExists(tixiHandle, xpath + "/structuralProfileUID");
            tixi::TixiSaveElement(tixiHandle, xpath + "/structuralProfileUID", *m_structuralProfileUID_choice2_2);
        }
        else {
            if (tixi::TixiCheckElement(tixiHandle, xpath + "/structuralProfileUID")) {
                tixi::TixiRemoveElement(tixiHandle, xpath + "/structuralProfileUID");
            }
        }

        // write element pointProperties
        tixi::TixiSaveElements(tixiHandle, xpath + "/pointProperties", m_pointProperties_choice2_2);

        // write element referencePointUID
        if (m_referencePointUID_choice2_2) {
            tixi::TixiCreateElementIfNotExists(tixiHandle, xpath + "/referencePointUID");
            tixi::TixiSaveElement(tixiHandle, xpath + "/referencePointUID", *m_referencePointUID_choice2_2);
        }
        else {
            if (tixi::TixiCheckElement(tixiHandle, xpath + "/referencePointUID")) {
                tixi::TixiRemoveElement(tixiHandle, xpath + "/referencePointUID");
            }
        }

    }

    bool CPACSProfileBasedStructuralElement::ValidateChoices() const
    {
        return
        (
            (
                (
                    // mandatory elements of this choice must be there
                    true // m_globalBeamProperties_choice1 is optional in choice
                    &&
                    // elements of other choices must not be there
                    !(
                        !m_sheetProperties_choice2.empty()
                        ||
                        m_standardProfileType_choice2_1.is_initialized()
                        ||
                        m_structuralProfileUID_choice2_2.is_initialized()
                        ||
                        !m_pointProperties_choice2_2.empty()
                        ||
                        m_referencePointUID_choice2_2.is_initialized()
                    )
                )
                +
                (
                    // mandatory elements of this choice must be there
                    !m_sheetProperties_choice2.empty()
                    &&
                    (
                        (
                            // mandatory elements of this choice must be there
                            m_standardProfileType_choice2_1.is_initialized()
                            &&
                            // elements of other choices must not be there
                            !(
                                m_structuralProfileUID_choice2_2.is_initialized()
                                ||
                                !m_pointProperties_choice2_2.empty()
                                ||
                                m_referencePointUID_choice2_2.is_initialized()
                            )
                        )
                        +
                        (
                            // mandatory elements of this choice must be there
                            m_structuralProfileUID_choice2_2.is_initialized()
                            &&
                            true // m_pointProperties_choice2_2 is optional in choice
                            &&
                            true // m_referencePointUID_choice2_2 is optional in choice
                            &&
                            // elements of other choices must not be there
                            !(
                                m_standardProfileType_choice2_1.is_initialized()
                            )
                        )
                        == 1
                    )
                    &&
                    // elements of other choices must not be there
                    !(
                        m_globalBeamProperties_choice1.is_initialized()
                    )
                )
                == 1
            )
        )
        ;
    }

    const std::string& CPACSProfileBasedStructuralElement::GetUID() const
    {
        return m_uID;
    }

    void CPACSProfileBasedStructuralElement::SetUID(const std::string& value)
    {
        if (m_uidMgr) {
            m_uidMgr->TryUnregisterObject(m_uID);
            m_uidMgr->RegisterObject(value, *this);
        }
        m_uID = value;
    }

    const boost::optional<std::string>& CPACSProfileBasedStructuralElement::GetName() const
    {
        return m_name;
    }

    void CPACSProfileBasedStructuralElement::SetName(const boost::optional<std::string>& value)
    {
        m_name = value;
    }

    const boost::optional<std::string>& CPACSProfileBasedStructuralElement::GetDescription() const
    {
        return m_description;
    }

    void CPACSProfileBasedStructuralElement::SetDescription(const boost::optional<std::string>& value)
    {
        m_description = value;
    }

    const boost::optional<CPACSTransformation2D>& CPACSProfileBasedStructuralElement::GetTransformation() const
    {
        return m_transformation;
    }

    boost::optional<CPACSTransformation2D>& CPACSProfileBasedStructuralElement::GetTransformation()
    {
        return m_transformation;
    }

    const boost::optional<CPACSGlobalBeamProperties>& CPACSProfileBasedStructuralElement::GetGlobalBeamProperties_choice1() const
    {
        return m_globalBeamProperties_choice1;
    }

    boost::optional<CPACSGlobalBeamProperties>& CPACSProfileBasedStructuralElement::GetGlobalBeamProperties_choice1()
    {
        return m_globalBeamProperties_choice1;
    }

    const std::vector<unique_ptr<CPACSMaterialDefinitionForProfileBased> >& CPACSProfileBasedStructuralElement::GetSheetProperties_choice2() const
    {
        return m_sheetProperties_choice2;
    }

    std::vector<unique_ptr<CPACSMaterialDefinitionForProfileBased> >& CPACSProfileBasedStructuralElement::GetSheetProperties_choice2()
    {
        return m_sheetProperties_choice2;
    }

    const boost::optional<CPACSProfileBasedStructuralElement_standardProfileType>& CPACSProfileBasedStructuralElement::GetStandardProfileType_choice2_1() const
    {
        return m_standardProfileType_choice2_1;
    }

    void CPACSProfileBasedStructuralElement::SetStandardProfileType_choice2_1(const boost::optional<CPACSProfileBasedStructuralElement_standardProfileType>& value)
    {
        m_standardProfileType_choice2_1 = value;
    }

    const boost::optional<std::string>& CPACSProfileBasedStructuralElement::GetStructuralProfileUID_choice2_2() const
    {
        return m_structuralProfileUID_choice2_2;
    }

    void CPACSProfileBasedStructuralElement::SetStructuralProfileUID_choice2_2(const boost::optional<std::string>& value)
    {
        m_structuralProfileUID_choice2_2 = value;
    }

    const std::vector<unique_ptr<CPACSMaterialDefinitionForProfileBasedPoint> >& CPACSProfileBasedStructuralElement::GetPointProperties_choice2_2() const
    {
        return m_pointProperties_choice2_2;
    }

    std::vector<unique_ptr<CPACSMaterialDefinitionForProfileBasedPoint> >& CPACSProfileBasedStructuralElement::GetPointProperties_choice2_2()
    {
        return m_pointProperties_choice2_2;
    }

    const boost::optional<std::string>& CPACSProfileBasedStructuralElement::GetReferencePointUID_choice2_2() const
    {
        return m_referencePointUID_choice2_2;
    }

    void CPACSProfileBasedStructuralElement::SetReferencePointUID_choice2_2(const boost::optional<std::string>& value)
    {
        m_referencePointUID_choice2_2 = value;
    }

    CPACSTransformation2D& CPACSProfileBasedStructuralElement::GetTransformation(CreateIfNotExistsTag)
    {
        if (!m_transformation)
            m_transformation = boost::in_place(m_uidMgr);
        return *m_transformation;
    }

    void CPACSProfileBasedStructuralElement::RemoveTransformation()
    {
        m_transformation = boost::none;
    }

    CPACSGlobalBeamProperties& CPACSProfileBasedStructuralElement::GetGlobalBeamProperties_choice1(CreateIfNotExistsTag)
    {
        if (!m_globalBeamProperties_choice1)
            m_globalBeamProperties_choice1 = boost::in_place(m_uidMgr);
        return *m_globalBeamProperties_choice1;
    }

    void CPACSProfileBasedStructuralElement::RemoveGlobalBeamProperties_choice1()
    {
        m_globalBeamProperties_choice1 = boost::none;
    }

    CPACSMaterialDefinitionForProfileBased& CPACSProfileBasedStructuralElement::AddSheetProperties_choice2()
    {
        m_sheetProperties_choice2.push_back(make_unique<CPACSMaterialDefinitionForProfileBased>());
        return *m_sheetProperties_choice2.back();
    }

    void CPACSProfileBasedStructuralElement::RemoveSheetProperties_choice2(CPACSMaterialDefinitionForProfileBased& ref)
    {
        for (std::size_t i = 0; i < m_sheetProperties_choice2.size(); i++) {
            if (m_sheetProperties_choice2[i].get() == &ref) {
                m_sheetProperties_choice2.erase(m_sheetProperties_choice2.begin() + i);
                return;
            }
        }
        throw CTiglError("Element not found");
    }

    CPACSMaterialDefinitionForProfileBasedPoint& CPACSProfileBasedStructuralElement::AddPointProperties_choice2_2()
    {
        m_pointProperties_choice2_2.push_back(make_unique<CPACSMaterialDefinitionForProfileBasedPoint>());
        return *m_pointProperties_choice2_2.back();
    }

    void CPACSProfileBasedStructuralElement::RemovePointProperties_choice2_2(CPACSMaterialDefinitionForProfileBasedPoint& ref)
    {
        for (std::size_t i = 0; i < m_pointProperties_choice2_2.size(); i++) {
            if (m_pointProperties_choice2_2[i].get() == &ref) {
                m_pointProperties_choice2_2.erase(m_pointProperties_choice2_2.begin() + i);
                return;
            }
        }
        throw CTiglError("Element not found");
    }

} // namespace generated
} // namespace tigl
