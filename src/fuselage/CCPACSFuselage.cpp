/* 
* Copyright (C) 2007-2013 German Aerospace Center (DLR/SC)
*
* Created: 2010-08-13 Markus Litz <Markus.Litz@dlr.de>
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
/**
* @file
* @brief  Implementation of CPACS fuselage handling routines.
*/
#include <cmath>
#include <iostream>

#include "tigl_config.h"

#include "CCPACSFuselage.h"
#include "CCPACSFuselageSegment.h"
#include "CCPACSFuselageStringerFramePosition.h"
#include "CCPACSConfiguration.h"
#include "CCPACSWingSegment.h"
#include "tiglcommonfunctions.h"
#include "CNamedShape.h"
#include "Debugging.h"
#include "CTiglCurveConnector.h"
#include "CTiglMakeLoft.h"
#include "CTiglBSplineAlgorithms.h"
#include "CCPACSFuselageSectionElement.h"
#include "CCPACSFuselageSectionElements.h"
#include "CCPACSFuselageSection.h"

#include "BRepOffsetAPI_ThruSections.hxx"
#include "BRepAlgoAPI_Fuse.hxx"
#include "ShapeFix_Shape.hxx"
#include "GProp_GProps.hxx"
#include "BRep_Tool.hxx"
#include "BRepTools.hxx"
#include "BRepGProp.hxx"
#include "BRepBuilderAPI_Transform.hxx"
#include "Geom_TrimmedCurve.hxx"
#include "TopoDS_Edge.hxx"
#include "BRepBuilderAPI_MakeEdge.hxx"
#include "BRepBuilderAPI_MakeFace.hxx"
#include "BRepBuilderAPI_MakeWire.hxx"
#include "GC_MakeSegment.hxx"
#include "BRepExtrema_DistShapeShape.hxx"
#include "ShapeFix_Wire.hxx"
#include "TopExp.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include <TopExp_Explorer.hxx>
#include <IntCurvesFace_Intersector.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepProj_Projection.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>

#include <math.h>       /* atan2 */


namespace tigl
{

CCPACSFuselage::CCPACSFuselage(CCPACSFuselages* parent, CTiglUIDManager* uidMgr)
    : generated::CPACSFuselage(parent, uidMgr)
    , CTiglRelativelyPositionedComponent(&m_parentUID, &m_transformation, &m_symmetry)
    , guideCurves(*this, &CCPACSFuselage::BuildGuideCurves) {
    Cleanup();
    if (parent->IsParent<CCPACSAircraftModel>())
        configuration = &parent->GetParent<CCPACSAircraftModel>()->GetConfiguration();
    else
        configuration = &parent->GetParent<CCPACSRotorcraftModel>()->GetConfiguration();
}

// Destructor
CCPACSFuselage::~CCPACSFuselage()
{
    Cleanup();
}

// Invalidates internal state
void CCPACSFuselage::Invalidate()
{
    loft.clear();
    guideCurves.clear();
    m_segments.Invalidate();
    if (m_positionings)
        m_positionings->Invalidate();
    if (m_structure)
        m_structure->Invalidate();
}

// Cleanup routine
void CCPACSFuselage::Cleanup()
{
    m_name = "";
    m_transformation.reset();

    // Calls ITiglGeometricComponent interface Reset to delete e.g. all childs.
    Reset();

    Invalidate();
}

// Read CPACS fuselage element
void CCPACSFuselage::ReadCPACS(const TixiDocumentHandle& tixiHandle, const std::string& fuselageXPath)
{
    Cleanup();

    generated::CPACSFuselage::ReadCPACS(tixiHandle, fuselageXPath);

    ConnectGuideCurveSegments();
}

// Returns the parent configuration
CCPACSConfiguration& CCPACSFuselage::GetConfiguration() const
{
    return *configuration;
}

std::string CCPACSFuselage::GetDefaultedUID() const
{
    return generated::CPACSFuselage::GetUID();
}

PNamedShape CCPACSFuselage::GetLoft(TiglCoordinateSystem cs) const
{
    PNamedShape loft = CTiglRelativelyPositionedComponent::GetLoft();
    if (!loft) {
        return loft;
    }

    if (cs == GLOBAL_COORDINATE_SYSTEM) {
        return loft;
    }
    else {
        // we want to modify the shape. we have to create a copy first
        loft = loft->DeepCopy();
        TopoDS_Shape transformedLoft = GetTransformationMatrix().Inverted().Transform(loft->Shape());
        loft->SetShape(transformedLoft);
        return loft;
    }
}

// Get section count
int CCPACSFuselage::GetSectionCount() const
{
    return m_sections.GetSectionCount();
}

// Returns the section for a given index
CCPACSFuselageSection& CCPACSFuselage::GetSection(int index) const
{
    return m_sections.GetSection(index);
}

// Get segment count
int CCPACSFuselage::GetSegmentCount() const
{
    return m_segments.GetSegmentCount();
}

// Returns the segment for a given index
CCPACSFuselageSegment& CCPACSFuselage::GetSegment(const int index)
{
    return m_segments.GetSegment(index);
}

const CCPACSFuselageSegment& CCPACSFuselage::GetSegment(const int index) const
{
    return m_segments.GetSegment(index);
}

// Returns the segment for a given uid
CCPACSFuselageSegment& CCPACSFuselage::GetSegment(std::string uid)
{
    return m_segments.GetSegment(uid);
}

// get short name for loft
std::string CCPACSFuselage::GetShortShapeName () const
{
    unsigned int findex = 0;
    for (int i = 1; i <= GetConfiguration().GetFuselageCount(); ++i) {
        tigl::CCPACSFuselage& f = GetConfiguration().GetFuselage(i);
        if (GetUID() == f.GetUID()) {
            findex = i;
            std::stringstream shortName;
            shortName << "F" << findex;
            return shortName.str();
        }
    }
    return "UNKNOWN";
}

void CCPACSFuselage::SetFaceTraits (PNamedShape loft) const
{
    int nFaces = GetNumberOfFaces(loft->Shape());
    bool hasSymmetryPlane = GetNumberOfEdges(m_segments.GetSegment(1).GetEndWire()) > 1;

    std::vector<std::string> names;
    names.push_back(loft->Name());
    names.push_back("symmetry");
    names.push_back("Front");
    names.push_back("Rear");

    // if we have a smooth surface, the whole fuslage is treatet as one segment
    int nSegments = this->GetSegmentCount();

    // the number of faces per segment depends on the number of guide curves and the existence of the symmetry plane
    int facesPerSegment = GetSegment(1).GetNumberOfLoftFaces();
    int remainingFaces = nFaces - facesPerSegment * nSegments;
    if (facesPerSegment == 0 || remainingFaces < 0 || remainingFaces > 2) {
        LOG(WARNING) << "Fuselage faces cannot be names properly";
        return;
    }

    int iFaceTotal = 0;
    int nSymmetryFaces = (int) hasSymmetryPlane;
    for (int iSegment = 0; iSegment < nSegments; ++iSegment) {
        for (int iFace = 0; iFace < facesPerSegment - nSymmetryFaces; ++iFace) {
            loft->FaceTraits(iFaceTotal++).SetName(names[0].c_str());
        }
        for (int iFace = 0; iFace < nSymmetryFaces; ++iFace) {
            loft->FaceTraits(iFaceTotal++).SetName(names[1].c_str());
        }
    }

    // set the caps
    int iFace = 2;
    for (;iFaceTotal < nFaces; ++iFaceTotal) {
        loft->FaceTraits(iFaceTotal).SetName(names[iFace++].c_str());
    }
}

// Builds a fused shape of all fuselage segments
PNamedShape CCPACSFuselage::BuildLoft() const
{
    TiglContinuity cont = m_segments.GetSegment(1).GetContinuity();
    Standard_Boolean smooth = (cont == ::C0? false : true);

    CTiglMakeLoft lofter;
    // add profiles
    for (int i=1; i <= m_segments.GetSegmentCount(); i++) {
        lofter.addProfiles(m_segments.GetSegment(i).GetStartWire());
    }
    lofter.addProfiles(m_segments.GetSegment(m_segments.GetSegmentCount()).GetEndWire());

    // add guides
    lofter.addGuides(GetGuideCurveWires());

    lofter.setMakeSolid(true);
    lofter.setMakeSmooth(smooth);

    TopoDS_Shape loftShape =  lofter.Shape();

    std::string loftName = GetUID();
    std::string loftShortName = GetShortShapeName();
    PNamedShape loft(new CNamedShape(loftShape, loftName.c_str(), loftShortName.c_str()));
    SetFaceTraits(loft);

    return loft;
}

// Get the positioning transformation for a given section index
boost::optional<CTiglTransformation> CCPACSFuselage::GetPositioningTransformation(const std::string &sectionUID)
{
    boost::optional<CTiglTransformation> ret;
    if (m_positionings)
        ret = m_positionings->GetPositioningTransformation(sectionUID);
    return ret;
}

// Gets a point on the given fuselage segment in dependence of a parameters eta and zeta with
// 0.0 <= eta <= 1.0 and 0.0 <= zeta <= 1.0. For eta = 0.0 the point lies on the start
// profile of the segment, for eta = 1.0 on the end profile of the segment. For zeta = 0.0
// the point is the start point of the profile wire, for zeta = 1.0 the last profile wire point.
gp_Pnt CCPACSFuselage::GetPoint(int segmentIndex, double eta, double zeta)
{
    return ((CCPACSFuselageSegment &) GetSegment(segmentIndex)).GetPoint(eta, zeta);
}


// Returns the volume of this fuselage
double CCPACSFuselage::GetVolume()
{
    const TopoDS_Shape fusedSegments = GetLoft()->Shape();

    // Calculate volume
    GProp_GProps System;
    BRepGProp::VolumeProperties(fusedSegments, System);
    double myVolume = System.Mass();
    return myVolume;
}

// Returns the circumference of the segment "segmentIndex" at a given eta
double CCPACSFuselage::GetCircumference(const int segmentIndex, const double eta)
{
    return static_cast<CCPACSFuselageSegment&>(GetSegment(segmentIndex)).GetCircumference(eta);
}
    
// Returns the surface area of this fuselage
double CCPACSFuselage::GetSurfaceArea()
{
    const PNamedShape& fusedSegments = GetLoft();

    // loop over all faces that are not symmetry, front or rear
    double myArea = 0.;

    TopTools_IndexedMapOfShape shapeMap;
    TopExp::MapShapes(fusedSegments->Shape(), TopAbs_FACE, shapeMap);
    for (int i = 1; i <= shapeMap.Extent(); ++i) {
        if (GetUID() == fusedSegments->GetFaceTraits(i-1).Name()) {
            const TopoDS_Face& f = TopoDS::Face(shapeMap(i));
            GProp_GProps System;
            BRepGProp::SurfaceProperties(f, System);
            myArea += System.Mass();
        }
    }

    // Calculate surface area
    return myArea;
}

// Returns the point where the distance between the selected fuselage and the ground is at minimum.
// The Fuselage could be turned with a given angle at at given axis, specified by a point and a direction.
gp_Pnt CCPACSFuselage::GetMinumumDistanceToGround(gp_Ax1 RAxis, double angle)
{

    TopoDS_Shape fusedFuselage = GetLoft()->Shape();

    // now rotate the fuselage
    gp_Trsf myTrsf;
    myTrsf.SetRotation(RAxis, angle * M_PI / 180.);
    BRepBuilderAPI_Transform xform(fusedFuselage, myTrsf);
    fusedFuselage = xform.Shape();

    // build cutting plane for intersection
    // We move the "ground" to "-1000" to be sure it is _under_ the fuselage
    gp_Pnt p1(-1.0e7, -1.0e7, -1000);
    gp_Pnt p2( 1.0e7, -1.0e7, -1000);
    gp_Pnt p3( 1.0e7,  1.0e7, -1000);
    gp_Pnt p4(-1.0e7,  1.0e7, -1000);

    Handle(Geom_TrimmedCurve) shaft_line1 = GC_MakeSegment(p1,p2);
    Handle(Geom_TrimmedCurve) shaft_line2 = GC_MakeSegment(p2,p3);
    Handle(Geom_TrimmedCurve) shaft_line3 = GC_MakeSegment(p3,p4);
    Handle(Geom_TrimmedCurve) shaft_line4 = GC_MakeSegment(p4,p1);

    TopoDS_Edge shaft_edge1 = BRepBuilderAPI_MakeEdge(shaft_line1);
    TopoDS_Edge shaft_edge2 = BRepBuilderAPI_MakeEdge(shaft_line2);
    TopoDS_Edge shaft_edge3 = BRepBuilderAPI_MakeEdge(shaft_line3);
    TopoDS_Edge shaft_edge4 = BRepBuilderAPI_MakeEdge(shaft_line4);

    TopoDS_Wire shaft_wire = BRepBuilderAPI_MakeWire(shaft_edge1, shaft_edge2, shaft_edge3, shaft_edge4);
    TopoDS_Face shaft_face = BRepBuilderAPI_MakeFace(shaft_wire);

    // calculate extrema
    BRepExtrema_DistShapeShape extrema(fusedFuselage, shaft_face);
    extrema.Perform();

    return extrema.PointOnShape1(1);
}

// Get the guide curve with a given UID
CCPACSGuideCurve& CCPACSFuselage::GetGuideCurveSegment(std::string uid)
{
    return const_cast<CCPACSGuideCurve&>(static_cast<const CCPACSFuselage&>(*this).GetGuideCurveSegment(uid));
}

const CCPACSGuideCurve& CCPACSFuselage::GetGuideCurveSegment(std::string uid) const
{
    for (int i = 1; i <= m_segments.GetSegmentCount(); i++) {
        const CCPACSFuselageSegment& segment = m_segments.GetSegment(i);

        if (!segment.GetGuideCurves()) {
            continue;
        }

        if (segment.GetGuideCurves()->GuideCurveExists(uid)) {
            return segment.GetGuideCurves()->GetGuideCurve(uid);
        }
    }
    throw tigl::CTiglError("Guide Curve with UID " + uid + " does not exists", TIGL_ERROR);
}

const TopoDS_Compound &CCPACSFuselage::GetGuideCurveWires() const
{
    return *guideCurves;
}

std::vector<gp_Pnt> CCPACSFuselage::GetGuideCurvePoints() const
{
    std::vector<gp_Pnt> points;

    // connect the belonging guide curve segments
    for (int isegment = 1; isegment <= GetSegmentCount(); ++isegment) {
        const CCPACSFuselageSegment& segment = m_segments.GetSegment(isegment);

        if (!segment.GetGuideCurves()) {
            continue;
        }

        const CCPACSGuideCurves& segmentCurves = *segment.GetGuideCurves();
        for (int iguide = 1; iguide <=  segmentCurves.GetGuideCurveCount(); ++iguide) {
            const CCPACSGuideCurve& curve = segmentCurves.GetGuideCurve(iguide);
            std::vector<gp_Pnt> curPoints = curve.GetCurvePoints();
            points.insert(points.end(), curPoints.begin(), curPoints.end());
        }
    }
    return points;
}

gp_Lin CCPACSFuselage::Intersection(gp_Pnt pRef, double angleRef) const
{
    // to have a left-handed coordinates system for the intersection computation (see documentation)
    const gp_Ax1 xAxe(pRef, gp_Dir(1, 0, 0));
    const gp_Dir zReference(0, 0, 1);
    const gp_Dir angleDir = zReference.Rotated(xAxe, angleRef);

    // build a line to position the intersection with the fuselage shape
    gp_Lin line(pRef, angleDir);

    // fuselage loft
    const TopoDS_Shape loft = GetLoft(FUSELAGE_COORDINATE_SYSTEM)->Shape();

    // get the list of shape from the fuselage shape
    TopExp_Explorer exp;
    for (exp.Init(loft, TopAbs_FACE); exp.More(); exp.Next()) {
        IntCurvesFace_Intersector intersection(TopoDS::Face(exp.Current()), 0.1); // intersection builder
        intersection.Perform(line, 0, std::numeric_limits<Standard_Real>::max());
        if (intersection.IsDone() && intersection.NbPnt() > 0) {
            gp_Lin result(intersection.Pnt(1), line.Direction());
            // return the line with the point on the fuselage as the origin, and the previous line's direction
            return result;
        }
    }

    TRACE_POINT(debug);
    debug.dumpShape(loft, "loft");
    debug.dumpShape(BRepBuilderAPI_MakeEdge(pRef, pRef.XYZ() + angleDir.XYZ() * 1000), "line");

    throw std::logic_error("Error computing intersection line");
}

gp_Lin CCPACSFuselage::Intersection(const CCPACSFuselageStringerFramePosition& pos) const
{
    const gp_Pnt pRef        = pos.GetRefPoint();
    const double angleRefRad = (M_PI / 180.) * pos.GetReferenceAngle();
    return Intersection(pRef, angleRefRad);
}

namespace
{
    TopoDS_Wire project(TopoDS_Shape wireOrEdge, BRepProj_Projection& proj, DebugScope& debug)
    {
        BRepBuilderAPI_MakeWire wireBuilder;
        for (; proj.More(); proj.Next())
            wireBuilder.Add(proj.Current());

        TopTools_ListOfShape wireList;
        BuildWiresFromConnectedEdges(proj.Shape(), wireList);

        if (wireList.Extent() == 0) {
            debug.addShape(proj.Shape(), "projection");

            throw CTiglError("Projection returned no wires");
        }
        if (wireList.Extent() == 1)
            return TopoDS::Wire(wireList.First());

        // select the wire which is closest to the wire we projected
        for (TopTools_ListIteratorOfListOfShape it(wireList); it.More(); it.Next()) {
            const TopoDS_Wire w                = TopoDS::Wire(it.Value());
            const gp_Pnt wStart     = GetFirstPoint(w);
            const gp_Pnt wEnd       = GetLastPoint(w);
            const gp_Pnt inputStart = GetFirstPoint(wireOrEdge);
            const gp_Pnt inputEnd   = GetLastPoint(wireOrEdge);

            const double pointEqualEpsilon = 1e-4;
            if ((wStart.IsEqual(inputStart, pointEqualEpsilon) && wEnd.IsEqual(inputEnd, pointEqualEpsilon)) ||
                (wEnd.IsEqual(inputStart, pointEqualEpsilon) && wStart.IsEqual(inputEnd, pointEqualEpsilon))) {
                return w;
            }
        }

        TopoDS_Compound c;
        TopoDS_Builder b;
        b.MakeCompound(c);
        for (TopTools_ListIteratorOfListOfShape it(wireList); it.More(); it.Next()) {
            b.Add(c, it.Value());
        }
        debug.addShape(proj.Shape(), "projection");
        debug.addShape(c, "wireList");

        // give up
        throw CTiglError("Failed to project wire/edge onto fuselage");
    }
}

TopoDS_Wire CCPACSFuselage::projectConic(TopoDS_Shape wireOrEdge, gp_Pnt origin) const
{
    const TopoDS_Shape loft = GetLoft(FUSELAGE_COORDINATE_SYSTEM)->Shape();

    DEBUG_SCOPE(debug);
    debug.addShape(wireOrEdge, "wireOrEdge");
    debug.addShape(loft, "loft");
    debug.addShape(BRepBuilderAPI_MakeVertex(origin), "origin");

    BRepProj_Projection proj(wireOrEdge, loft, origin);
    return project(wireOrEdge, proj, debug);
}

TopoDS_Wire CCPACSFuselage::projectParallel(TopoDS_Shape wireOrEdge, gp_Dir direction) const
{
    const TopoDS_Shape loft = GetLoft(FUSELAGE_COORDINATE_SYSTEM)->Shape();

    const TopoDS_Shape directionLine = BRepBuilderAPI_MakeEdge(
        BRepBuilderAPI_MakeVertex(gp_Pnt(0, 0, 0)).Vertex(),
        BRepBuilderAPI_MakeVertex(gp_Pnt(direction.XYZ() * 1000)).Vertex()
    ).Shape();

    DEBUG_SCOPE(debug);
    debug.addShape(wireOrEdge, "wireOrEdge");
    debug.addShape(loft, "loft");
    debug.addShape(directionLine, "direction");

    BRepProj_Projection proj(wireOrEdge, loft, direction);
    return project(wireOrEdge, proj, debug);
}

void CCPACSFuselage::BuildGuideCurves(TopoDS_Compound& cache) const
{
    std::map<double, const CCPACSGuideCurve*> roots;

    // get section centers for the centripetal parametrization
    std::vector<gp_Pnt> sectionCenters(GetSegmentCount()+1);

    // get center of inner section of first segment
    const CCPACSFuselageSegment& innerSegment = m_segments.GetSegment(1);
    sectionCenters[0] = innerSegment.GetTransformedProfileOriginStart();
    
    // find roots and connect the belonging guide curve segments
    for (int isegment = 1; isegment <= GetSegmentCount(); ++isegment) {
        const CCPACSFuselageSegment& segment = m_segments.GetSegment(isegment);

        if (!segment.GetGuideCurves()) {
            continue;
        }

        // get center of outer section
        sectionCenters[isegment] = segment.GetTransformedProfileOriginEnd();

        const CCPACSGuideCurves& segmentCurves = *segment.GetGuideCurves();
        for (int iguide = 1; iguide <=  segmentCurves.GetGuideCurveCount(); ++iguide) {
            const CCPACSGuideCurve& curve = segmentCurves.GetGuideCurve(iguide);
            if (!curve.GetFromGuideCurveUID_choice1()) {
                // this is a root curve
                double relCirc= *curve.GetFromRelativeCircumference_choice2();
                //TODO: determine if half fuselage or not. If not
                //the guide curve at relCirc=1 should be inserted at relCirc=0
                roots.insert(std::make_pair(relCirc, &curve));
            }
        }
    }

    // get the parameters at the section centers
    std::vector<double> sectionParams = CTiglBSplineAlgorithms::computeParamsBSplineCurve(OccArray(sectionCenters), 0., 1., 0.5);

    // connect guide curve segments to a spline with given continuity conditions and tangents
    CTiglCurveConnector connector(roots, sectionParams);
    cache = connector.GetConnectedGuideCurves();
}

void CCPACSFuselage::ConnectGuideCurveSegments(void)
{
    for (int isegment = 1; isegment <= GetSegmentCount(); ++isegment) {
        CCPACSFuselageSegment& segment = GetSegment(isegment);

        if (!segment.GetGuideCurves()) {
            continue;
        }

        CCPACSGuideCurves& curves = *segment.GetGuideCurves();
        for (int icurve = 1; icurve <= curves.GetGuideCurveCount(); ++icurve) {
            CCPACSGuideCurve& curve = curves.GetGuideCurve(icurve);
            if (!curve.GetFromRelativeCircumference_choice2()) {
                std::string fromUID = *curve.GetFromGuideCurveUID_choice1();
                CCPACSGuideCurve& fromCurve = GetGuideCurveSegment(fromUID);
                curve.SetFromRelativeCircumference_choice2(fromCurve.GetToRelativeCircumference());
            }
        }
    }
}

TopoDS_Shape transformFuselageProfileGeometry(const CTiglTransformation& fuselTransform, const CTiglFuselageConnection& connection, const TopoDS_Shape& shape)
{
    // Do section element transformation on points
    tigl::CTiglTransformation trafo = connection.GetSectionElementTransformation();

    // Do section transformations
    trafo.PreMultiply(connection.GetSectionTransformation());

    // Do positioning transformations
    boost::optional<CTiglTransformation> posTrans = connection.GetPositioningTransformation();
    if (posTrans) {
        trafo.PreMultiply(*posTrans);
    }

    trafo.PreMultiply(fuselTransform);

    return trafo.Transform(shape);

}

double CCPACSFuselage::GetLength()
{
    std::string noise = GetNoiseUID();
    std::string tail  = GetTailUID();
    return GetLengthBetween(noise, tail);
}

double CCPACSFuselage::GetLengthBetween(const std::string& startElementUID, const std::string& endElementUID)
{
    std::map<std::string, CTiglPoint>& centers = GetGraph().getCenters(); // is it passs by reference? yes // we do not use GetElementsCenters for performance reason
    CTiglPoint delta                          = centers[endElementUID] - centers[startElementUID];
    return delta.norm2();
}

std::map<std::string, CTiglPoint> CCPACSFuselage::GetElementsCenters()
{

    std::map<std::string, CTiglPoint> centers; // center of the element

    CCPACSFuselageSegments& segments = GetSegments();

    std::string tempFromElmUID, tempToElmUID;
    CTiglPoint fromCenter, toCenter;
    TopoDS_Shape curve;
    gp_Pnt centerPoint;
    for (int i = 1; i <= segments.GetSegmentCount(); i++) {
        CCPACSFuselageSegment& seg = segments.GetSegment(i);
        tempFromElmUID             = seg.GetStartSectionElementUID();
        tempToElmUID               = seg.GetEndSectionElementUID();

        // get the center point of elements for later
        curve                   = seg.GetStartWire();
        centerPoint             = GetCenterOfMass(curve);
        fromCenter              = CTiglPoint(centerPoint.XYZ());
        centers[tempFromElmUID] = fromCenter; // overwrite if existing or create if not

        curve                 = seg.GetEndWire();
        centerPoint           = GetCenterOfMass(curve);
        toCenter              = CTiglPoint(centerPoint.XYZ());
        centers[tempToElmUID] = toCenter;
    }

    return centers;
}

std::vector<std::string> CCPACSFuselage::GetCreatorGraph()
{
    std::vector<std::string> r = GetGraph().getCreatorGraph();  // should copy
    return r;

}

std::string CCPACSFuselage::GetNoiseUID()
{
    return GetGraph().getNoiseUID();
}

std::string CCPACSFuselage::GetTailUID()
{
    return GetGraph().getTailUID();
}

void CCPACSFuselage::SetLength(double newLength){
    std::string noise = GetNoiseUID();
    std::string tail = GetTailUID();
    SetLengthBetween(noise, tail, newLength);
}



void CCPACSFuselage::SetLengthBetween(const std::string& startElement, const std::string& endElement,
                                        double newPartialLength)
{

    // todo verify input order
    std::vector<std::string>& tempGraph = GetGraph().getCreatorGraph();

    /*
         * Divide the elements in 3 categories:
         * 1) Elements before start that need not to be modified
         * 2) Elements between that need to create the partial length
         * 3) Elements after end that need to be shifted has the last between element
         */
    std::vector<std::string> elementsBetween = GetGraph().getElementsInBetween(startElement, endElement);
    std::vector<std::string> elementsBefore = GetGraph().getElementsBefore(startElement);
    std::vector<std::string> elementsAfter = GetGraph().getElementsAfter(endElement);

    if( elementsBetween.size() < 2 ){
        throw CTiglError("CCPACSFuselage::SetLengthBetween: At least two elements should be contains between stratUID and endUID");
    }



    /*
         * BETWEEN ELEMENT SCALING
         *
         * This part follow basically these steps:
         *
         * 1) Computation of the affine transformations needed to perform the desired effect. The desired effect can be perform as:
         *         a) Put the start on the origin
         *         b) Rotation to get the end on the X axis
         *         c) Perform a scaling on X to obtain the desired length value
         *         d) inverse of of b) to put the fuselage in the right direction
         *         e) inverse of a) to shift the fuselage to its origin place
         *
         *
         * 2) Compute the origin of each element and delta between the origin and the center point.
         *    This is done to cover the case of profiles that are shifted (origin of element != center point)
         *
         * 3) Compute the new center points and the new origin to get the wanted length
         *
         * 4) Find the new Transformation that elements should have to have these new origins and save in the memory
         *
         *
         */

    std::map<std::string, CTiglPoint> oldGlobalOrigin;
    std::map<std::string, CTiglPoint> newCenterPoints;
    std::map<std::string, CTiglPoint> newGlobalOrigin;

    // Get fuselage point in world coordinate
    std::map<std::string, CTiglPoint>& oldCenterPoints   = GetGraph().getCenters();
    CTiglPoint startP = oldCenterPoints[startElement];
    CTiglPoint endP   = oldCenterPoints[endElement];

    // bring StartP to Origin
    CTiglTransformation startToO;
    startToO.SetIdentity();
    startToO.AddTranslation(-startP.x, -startP.y, -startP.z);

    startP = startToO * startP;
    endP   = startToO * endP;


    // bring endP on the x axis
    // We perform a extrinsic rotation in the order Z Y X, so it should be equivalent to the intrinsic cpacs rotation
    // in the order X Y' Z''
    CTiglTransformation rotEndToX4d;
    rotEndToX4d.SetIdentity();
    double rotGradZ = atan2(endP.y, endP.x);
    double rotZ = CTiglTransformation::RadianToDegree(rotGradZ);
    rotEndToX4d.AddRotationZ(-rotZ);
    double rotGradY = atan2(endP.z, sqrt( (endP.x * endP.x)  + (endP.y * endP.y )));
    double rotY = CTiglTransformation::RadianToDegree(rotGradY);
    rotEndToX4d.AddRotationY(-rotY);

    endP = rotEndToX4d * endP;

    double oldPartialLength = GetLengthBetween(startElement, endElement);

    // Compute the needed scaling in x
    if (oldPartialLength == 0) {
        throw CTiglError("CCPACSFuselage::SetLengthBetween: the old length is 0, impossible to scale the length");
    }

    double xScale = newPartialLength / oldPartialLength;
    CTiglTransformation scaleM;
    scaleM.SetIdentity();
    scaleM.AddScaling(xScale, 1.0, 1.0);

    CTiglTransformation startToOI    = startToO.Inverted();
    CTiglTransformation rotEndToX4dI = rotEndToX4d.Inverted();

    // Get the origin of each element
    CTiglPoint origin(0, 0, 0);
    CTiglTransformation global;
    for (int i = 0; i < tempGraph.size(); i++) {
        global = GetGlobalTransformation(tempGraph[i]);
        oldGlobalOrigin[tempGraph[i]] = global * origin;
    }

    // Compute the new center point and the new origin of each element in Between
    CTiglTransformation totalTransformation = startToOI * rotEndToX4dI * scaleM * rotEndToX4d * startToO;
    CTiglPoint tempDelatOtoP;
    for (int i = 0; i < elementsBetween.size(); i++) {
        newCenterPoints[elementsBetween[i]] = totalTransformation * oldCenterPoints[elementsBetween[i]];
        tempDelatOtoP                       = oldCenterPoints[elementsBetween[i]] - oldGlobalOrigin[elementsBetween[i]];
        // delta between origin and the center point will not change because no scaling or rotation will be changed
        newGlobalOrigin[elementsBetween[i]] = newCenterPoints[elementsBetween[i]] - tempDelatOtoP;
    }

    // Compute the new transformation element of each element to be placed at the wanted orgin
    CTiglTransformation tempNewTransformationE;
    for (int i = 0; i < elementsBetween.size(); i++) {
        tempNewTransformationE =
            GetTransformToPlaceElementByTranslationAt(elementsBetween[i], newGlobalOrigin[elementsBetween[i]]);
        CCPACSFuselageSectionElement& element =
            GetUIDManager().ResolveObject<CCPACSFuselageSectionElement>(elementsBetween[i]);
        CCPACSTransformation& storedTransformation = element.GetTransformation();
        storedTransformation.setTransformationMatrix(tempNewTransformationE);
    }

    /*
     * SHIFT THE END OF THE FUSELAGE
    */

    CTiglPoint shiftEndElementG = newGlobalOrigin[endElement] - oldGlobalOrigin[endElement];
    for (int i = 0; i < elementsAfter.size(); i++) {

        newGlobalOrigin[elementsAfter[i]] = oldGlobalOrigin[elementsAfter[i]] + shiftEndElementG;
        tempNewTransformationE = GetTransformToPlaceElementByTranslationAt(elementsAfter[i], newGlobalOrigin[elementsAfter[i]]);
        CCPACSFuselageSectionElement& element = GetUIDManager().ResolveObject<CCPACSFuselageSectionElement>(elementsAfter[i]);
        CCPACSTransformation& storedTransformation = element.GetTransformation();
        storedTransformation.setTransformationMatrix(tempNewTransformationE);

    }


    // UPDATE THE FUSELAGE STRUCTURE

    GetGraph().clear(); // the structure has changed -> so we need to reconstruct the graph
    GetConfiguration().WriteCPACS(GetConfiguration().GetUID());
    Invalidate(); // to force the rebuild of the loft? yes
}

std::vector<CTiglTransformation> CCPACSFuselage::GetTransformationChain(const std::string& elementUID)
{
    std::vector<CTiglTransformation> result;
    CCPACSFuselageSectionElement& element = GetUIDManager().ResolveObject<CCPACSFuselageSectionElement>(elementUID);
    CTiglTransformation elementT          = element.GetSectionElementTransformation();
    CCPACSFuselageSection* section        = element.GetParent()->GetParent();
    CTiglTransformation sectionT          = section->GetSectionTransformation();
    CTiglTransformation positioningT;
    positioningT.SetIdentity();
    boost::optional<CTiglTransformation> positioningOp = GetPositioningTransformation(section->GetUID());
    if (positioningOp) {
        positioningT = positioningOp.get();
    }
    CTiglTransformation fuselageT = GetTransformationMatrix();
    result.push_back(elementT);
    result.push_back(sectionT);
    result.push_back(positioningT);
    result.push_back(fuselageT);

    return result;
}

CTiglTransformation CCPACSFuselage::GetGlobalTransformation(const std::string& elementUID)
{
    std::vector<CTiglTransformation> chain = GetTransformationChain(elementUID);
    CTiglTransformation global             = chain[3] * chain[2] * chain[1] * chain[0];
    return global;
}

CTiglTransformation CCPACSFuselage::GetTransformToPlaceElementByTranslationAt(const std::string& elementUID,
                                                                              const CTiglPoint& wantedOriginP)
{

    /* We search a new E' such that:
         * w = F*P*S*E'*0  (At start the origin is at 0,0,0)
         * S^-1 * P^-1 * F ^-1 w = E'0
         * A w = E'0
         * w' = E' * 0
         * w' = T' *R' * S' * 0 ( we decompose the E' into its scaling, rotation and translation, remark the S has not the same meaning as above)
         * w' = T' * 0 ( Because scaling and rotation has no effect on 0)
         * so we can deduce the wanted E'
         * E' = T' * R * S where S and R are the original scaling and rotation of the original E.
         *
         */

    std::vector<CTiglTransformation> chain = GetTransformationChain(elementUID);

    CTiglTransformation a = chain[1].Inverted() * chain[2].Inverted() * chain[3].Inverted();
    CTiglPoint wp         = a * wantedOriginP;

    CTiglTransformation ep = chain[0];

    ep.SetValue(0, 3, 0);
    ep.SetValue(1, 3, 0);
    ep.SetValue(2, 3, 0);
    ep.SetValue(3, 3, 1);

    ep.AddTranslation(wp.x, wp.y, wp.z);

    // check if it is correct
    CTiglPoint o(0, 0, 0);
    CTiglPoint check = chain[3] * chain[2] * chain[1] * ep * o;

    if (check.distance2(wantedOriginP) > 0.001) {
        throw CTiglError("CCPACSFuselage::GetTransformToPlaceElementByTranslationAt: Something go wrong!");
    }
    return ep;
}

std::map<std::string, double> CCPACSFuselage::GetCircumferenceOfElements()
{

    std::map<std::string, double> circumferences;

    CCPACSFuselageSegments& segments = GetSegments();

    std::string tempFromElmUID, tempToElmUID;
    CTiglPoint fromCenter, toCenter;
    TopoDS_Shape curve;
    gp_Pnt centerPoint;
    GProp_GProps System;
    for (int i = 1; i <= segments.GetSegmentCount(); i++) {
        CCPACSFuselageSegment& seg = segments.GetSegment(i);
        tempFromElmUID             = seg.GetStartSectionElementUID();
        tempToElmUID               = seg.GetEndSectionElementUID();

        // We do not want to use this function GetCircumferenece function because they are to slow!
        //circumferences[tempFromElmUID] = GetCircumference(i, 0);
        //circumferences[tempToElmUID]   = GetCircumference(i, 1);
        curve                   = seg.GetStartWire();
        BRepGProp::LinearProperties(curve,System);
        circumferences[tempFromElmUID] = System.Mass();

        curve                   = seg.GetEndWire();
        BRepGProp::LinearProperties(curve,System);
        circumferences[tempToElmUID] = System.Mass();
    }

    return circumferences;
}

double CCPACSFuselage::GetMaximalCircumferenceOfElements()
{
    std::string noise = GetNoiseUID();
    std::string tail = GetTailUID();
    return GetMaximalCircumferenceOfElementsBetween(noise, tail);
}


double CCPACSFuselage::GetMaximalCircumferenceOfElementsBetween(std::string startElementUID, std::string endElementUID){

    std::vector<std::string> elementsBetween = GetGraph().getElementsInBetween(startElementUID, endElementUID);

    if (elementsBetween.size() < 1) {
        throw CTiglError("CCPACSFuselage::GetMaximalCircumferenceOfElementsBetween: impossible to get the start and the end correctly");
    }

    std::map<std::string, double>& circumferences = GetGraph().getCircumferences();

    std::string maxUID      = "";
    double maxCircumference = std::numeric_limits<double>::min();

    for(int i = 0; i < elementsBetween.size(); i++ ){
        if( circumferences[elementsBetween[i]]> maxCircumference){
            maxUID = elementsBetween[i];
            maxCircumference = circumferences[elementsBetween[i]];
        }
    }

    return maxCircumference;
}



void CCPACSFuselage::ScaleCircumferenceOfElements(std::vector<std::string> elementsToScale, double scaleFactor)
{

    /*
         *  The idea is to bring the element to its center and scale it in world coordinate system.
         *  So we get:
         *  FPSE' = Ti*S*T*FPSE
         *  E' =  Si*Pi*Fi*Ti*S*T*FPSE
         */

    CTiglTransformation centerToOriginM, centerToOriginMI, scaleM, newE;

    scaleM.SetIdentity();
    scaleM.AddScaling(scaleFactor, scaleFactor, scaleFactor);

    std::map<std::string, CTiglPoint>& centers = GetGraph().getCenters();
    std::vector<CTiglTransformation> chain;
    std::string tempUID;
    for (int i = 0; i < elementsToScale.size(); i++) {
        tempUID = elementsToScale[i];
        chain   = GetTransformationChain(tempUID);
        centerToOriginM.SetIdentity();
        centerToOriginM.AddTranslation(-centers[tempUID].x, -centers[tempUID].y, -centers[tempUID].z);
        centerToOriginMI = centerToOriginM.Inverted();
        newE             = chain[1].Inverted() * chain[2].Inverted() * chain[3].Inverted() * centerToOriginMI * scaleM *
               centerToOriginM * chain[3] * chain[2] * chain[1] * chain[0];
        CCPACSFuselageSectionElement& element =
                GetUIDManager().ResolveObject<CCPACSFuselageSectionElement>(tempUID);
        CCPACSTransformation& storedTransformation = element.GetTransformation();
        storedTransformation.setTransformationMatrix(newE);
    }


    // Update the tigl object
    GetGraph().clear();
    GetConfiguration().WriteCPACS(GetConfiguration().GetUID());
    Invalidate();
}

void CCPACSFuselage::SetMaximalCircumferenceOfElementsBetween(std::string startUID, std::string endUID,
                                                            double newMaximalCircumference)
{

    if (newMaximalCircumference < 0) {
        throw CTiglError("setFuselageMaximalCircumference: the newMaximalCricumference should be bigger than 0");
    }

    double oldMaximalCircumference = GetMaximalCircumferenceOfElementsBetween(startUID, endUID);

    double scaleFactor = newMaximalCircumference / oldMaximalCircumference;

    std::vector<std::string> elementsBetween = GetGraph().getElementsInBetween(startUID, endUID);

    if (elementsBetween.size() < 1) {
        throw CTiglError("CCPACSFuselage::GetMaximalCircumferenceOfElementsBetween: impossible to get the start and the end correctly");
    }

    ScaleCircumferenceOfElements(elementsBetween, scaleFactor);
}

void CCPACSFuselage::SetMaximalCircumferenceOfElements(double newMaximalCircumference)
{
    std::string noiseUID = GetNoiseUID();
    std::string tailUID = GetTailUID();
    SetMaximalCircumferenceOfElementsBetween(noiseUID, tailUID, newMaximalCircumference);

}

FuselageGraph& CCPACSFuselage::GetGraph()
{
    if( ! graph.isBuild()){
       graph.build(this);
    }
    return graph;
}

} // end namespace tigl
