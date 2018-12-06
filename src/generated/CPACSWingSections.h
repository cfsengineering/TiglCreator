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

#pragma once

#include <string>
#include <tixi.h>
#include <typeinfo>
#include <vector>
#include "CTiglError.h"
#include "tigl_internal.h"
#include "UniquePtr.h"

namespace tigl
{
class CTiglUIDManager;
class CCPACSWingSection;
class CCPACSEnginePylon;
class CCPACSWing;

namespace generated
{
    // This class is used in:
    // CPACSEnginePylon
    // CPACSWing

    // generated from /xsd:schema/xsd:complexType[951]
    /// @brief Sections of the wing.
    /// 
    /// WingSections type, containing all the sections of the
    /// wing.
    /// 
    class CPACSWingSections
    {
    public:
        TIGL_EXPORT CPACSWingSections(CCPACSEnginePylon* parent, CTiglUIDManager* uidMgr);
        TIGL_EXPORT CPACSWingSections(CCPACSWing* parent, CTiglUIDManager* uidMgr);

        TIGL_EXPORT virtual ~CPACSWingSections();

        template<typename P>
        bool IsParent() const
        {
            return m_parentType != NULL && *m_parentType == typeid(P);
        }

        template<typename P>
        P* GetParent()
        {
#ifdef HAVE_STDIS_SAME
            static_assert(std::is_same<P, CCPACSEnginePylon>::value || std::is_same<P, CCPACSWing>::value, "template argument for P is not a parent class of CPACSWingSections");
#endif
            if (!IsParent<P>()) {
                throw CTiglError("bad parent");
            }
            return static_cast<P*>(m_parent);
        }

        template<typename P>
        const P* GetParent() const
        {
#ifdef HAVE_STDIS_SAME
            static_assert(std::is_same<P, CCPACSEnginePylon>::value || std::is_same<P, CCPACSWing>::value, "template argument for P is not a parent class of CPACSWingSections");
#endif
            if (!IsParent<P>()) {
                throw CTiglError("bad parent");
            }
            return static_cast<P*>(m_parent);
        }

        TIGL_EXPORT CTiglUIDManager& GetUIDManager();
        TIGL_EXPORT const CTiglUIDManager& GetUIDManager() const;

        TIGL_EXPORT virtual void ReadCPACS(const TixiDocumentHandle& tixiHandle, const std::string& xpath);
        TIGL_EXPORT virtual void WriteCPACS(const TixiDocumentHandle& tixiHandle, const std::string& xpath) const;

        TIGL_EXPORT virtual const std::vector<unique_ptr<CCPACSWingSection> >& GetSections() const;
        TIGL_EXPORT virtual std::vector<unique_ptr<CCPACSWingSection> >& GetSections();

        TIGL_EXPORT virtual CCPACSWingSection& AddSection();
        TIGL_EXPORT virtual void RemoveSection(CCPACSWingSection& ref);

    protected:
        void* m_parent;
        const std::type_info* m_parentType;

        CTiglUIDManager* m_uidMgr;

        std::vector<unique_ptr<CCPACSWingSection> > m_sections;

    private:
#ifdef HAVE_CPP11
        CPACSWingSections(const CPACSWingSections&) = delete;
        CPACSWingSections& operator=(const CPACSWingSections&) = delete;

        CPACSWingSections(CPACSWingSections&&) = delete;
        CPACSWingSections& operator=(CPACSWingSections&&) = delete;
#else
        CPACSWingSections(const CPACSWingSections&);
        CPACSWingSections& operator=(const CPACSWingSections&);
#endif
    };
} // namespace generated

// CPACSWingSections is customized, use type CCPACSWingSections directly
} // namespace tigl
