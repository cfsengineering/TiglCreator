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

#include <boost/optional.hpp>
#include <boost/utility/in_place_factory.hpp>
#include <CCPACSEtaIsoLine.h>
#include <string>
#include <tixi.h>
#include "CPACSCap.h"
#include "CPACSWeb.h"
#include "CreateIfNotExists.h"
#include "tigl_internal.h"

namespace tigl
{
class CTiglUIDManager;

namespace generated
{
    // This class is used in:
    // CPACSSparCells

    // generated from /xsd:schema/xsd:complexType[810]
    /// @brief Spar cell of the spar.
    /// 
    /// Within spar cells a special area of the spar is
    /// defined where different cross section and material properties
    /// shall be defined.
    /// The area of the spar is defined by using the
    /// parameters 'fromEta' and 'toEta'. The definition of the caps,
    /// webs and rotation is equivalent to the cross section definition
    /// of the complete spar.
    /// 
    class CPACSSparCell
    {
    public:
        TIGL_EXPORT CPACSSparCell(CTiglUIDManager* uidMgr);
        TIGL_EXPORT virtual ~CPACSSparCell();

        TIGL_EXPORT CTiglUIDManager& GetUIDManager();
        TIGL_EXPORT const CTiglUIDManager& GetUIDManager() const;

        TIGL_EXPORT virtual void ReadCPACS(const TixiDocumentHandle& tixiHandle, const std::string& xpath);
        TIGL_EXPORT virtual void WriteCPACS(const TixiDocumentHandle& tixiHandle, const std::string& xpath) const;

        TIGL_EXPORT virtual const std::string& GetUID() const;
        TIGL_EXPORT virtual void SetUID(const std::string& value);

        TIGL_EXPORT virtual const CCPACSEtaIsoLine& GetFromEta() const;
        TIGL_EXPORT virtual CCPACSEtaIsoLine& GetFromEta();

        TIGL_EXPORT virtual const CCPACSEtaIsoLine& GetToEta() const;
        TIGL_EXPORT virtual CCPACSEtaIsoLine& GetToEta();

        TIGL_EXPORT virtual const CPACSCap& GetUpperCap() const;
        TIGL_EXPORT virtual CPACSCap& GetUpperCap();

        TIGL_EXPORT virtual const CPACSCap& GetLowerCap() const;
        TIGL_EXPORT virtual CPACSCap& GetLowerCap();

        TIGL_EXPORT virtual const CPACSWeb& GetWeb1() const;
        TIGL_EXPORT virtual CPACSWeb& GetWeb1();

        TIGL_EXPORT virtual const boost::optional<CPACSWeb>& GetWeb2() const;
        TIGL_EXPORT virtual boost::optional<CPACSWeb>& GetWeb2();

        TIGL_EXPORT virtual const double& GetRotation() const;
        TIGL_EXPORT virtual void SetRotation(const double& value);

        TIGL_EXPORT virtual CPACSWeb& GetWeb2(CreateIfNotExistsTag);
        TIGL_EXPORT virtual void RemoveWeb2();

    protected:
        CTiglUIDManager* m_uidMgr;

        std::string               m_uID;

        /// Beginning (= inner border) of the spar cell.
        CCPACSEtaIsoLine          m_fromEta;

        /// Ending (= outer border) of the spar cell.
        CCPACSEtaIsoLine          m_toEta;

        CPACSCap                  m_upperCap;

        CPACSCap                  m_lowerCap;

        CPACSWeb                  m_web1;

        boost::optional<CPACSWeb> m_web2;

        /// The angle between the wing middle plane and
        /// web1. Default is 90 degrees. Positive rotation is around the
        /// spar axis heading along with the positive eta-axis.
        double                    m_rotation;

    private:
#ifdef HAVE_CPP11
        CPACSSparCell(const CPACSSparCell&) = delete;
        CPACSSparCell& operator=(const CPACSSparCell&) = delete;

        CPACSSparCell(CPACSSparCell&&) = delete;
        CPACSSparCell& operator=(CPACSSparCell&&) = delete;
#else
        CPACSSparCell(const CPACSSparCell&);
        CPACSSparCell& operator=(const CPACSSparCell&);
#endif
    };
} // namespace generated

// Aliases in tigl namespace
#ifdef HAVE_CPP11
using CCPACSSparCell = generated::CPACSSparCell;
#else
typedef generated::CPACSSparCell CCPACSSparCell;
#endif
} // namespace tigl
