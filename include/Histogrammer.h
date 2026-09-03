/**
 * @file Histogrammer.h
 * @brief Manager for creating and persisting multidimensional kinematic histograms.
 * @details Wraps RDataFrame's histogramming capabilities to support automatic 
 *          splitting, masking, and unfolding of N-D histograms into 1D slices.
 */

#pragma once

#include "KinematicsProcessor.h"
#include "PhysicsSelection.h"
#include "ConfigReaction.h"
#include "THnCombi.h" 

#include <THnSparse.h>
#include <TFile.h>
#include <TDirectory.h>
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <stdexcept>

namespace rad {
namespace histo {

    /**
     * @struct HistoDef
     * @brief Configuration storage for lazy 1D histogram initialization.
     */
    struct HistoDef {
        std::string name;
        std::string title;
        int nbins;
        double min;
        double max;
        std::string varBaseName;
    };

    /**
     * @struct HistoDef2D
     * @brief Configuration storage for lazy 2D histogram initialization.
     */
    struct HistoDef2D {
        std::string name;
        std::string title;
        int nbinsX;
        double xmin;
        double xmax;
        int nbinsY;
        double ymin;
        double ymax;
        std::string XvarBaseName;
        std::string YvarBaseName;
    };

    /**
     * @struct SplitDef
     * @brief Configuration for a single splitting axis (category) like Sector.
     */
    struct SplitDef {
        std::string name;    
        std::string baseCol; 
        int nbins;
        double min;
        double max;
    };

    /**
     * @class Histogrammer
     * @brief Coordinates N-Dimensional THnSparse generation and automatic unfolding.
     */
    class Histogrammer {
    public:
        /**
         * @brief Constructor binding to a specific Processor stream.
         * @param proc The processor (Rec, Sim, etc.) used to resolve variable names.
         * @param sel Optional pointer to PhysicsSelection for automatic masking.
         */
        Histogrammer(KinematicsProcessor& proc, PhysicsSelection* sel = nullptr);

        /**
         * @brief Defines a split axis (category) for subsequent histograms.
         * @param name Axis title (e.g. "Sector").
         * @param baseCol Variable name in RDF (e.g. "loc_sector").
         * @param nbins Number of bins.
         * @param min Axis min.
         * @param max Axis max.
         */
        void AddSplit(const std::string& name, const std::string& baseCol, int nbins, double min, double max);

        /**
         * @brief Queues a 1D histogram for booking.
         * @param name Output histogram base name.
         * @param title Histogram title.
         * @param nbins Number of bins for the main variable.
         * @param min Minimum value for the main variable.
         * @param max Maximum value for the main variable.
         * @param varBaseName The RDF column name (e.g. "mm2").
         */
        void Create(const std::string& name, const std::string& title, 
                    int nbins, double min, double max, 
                    const std::string& varBaseName);

        /**
         * @brief Queues a 2D histogram for booking.
         * @param name Output histogram base name.
         * @param title Histogram title.
         * @param nbinsX Number of X bins.
         * @param xmin X minimum.
         * @param xmax X maximum.
         * @param nbinsY Number of Y bins.
         * @param ymin Y minimum.
         * @param ymax Y maximum.
         * @param XvarBaseName The RDF column name for the X axis.
         * @param YvarBaseName The RDF column name for the Y axis.
         */
        void Create2D(const std::string& name, const std::string& title,
                      int nbinsX, double xmin, double xmax,
                      int nbinsY, double ymin, double ymax,
                      const std::string& XvarBaseName,
                      const std::string& YvarBaseName);
    
        /**
         * @brief Compiles definitions and Books actions in RDataFrame.
         */
        void Init();

        /**
         * @brief Writes and unfolds all histograms to a ROOT file.
         * @param filename Output filename.
         * @param option ROOT File option ("RECREATE" or "UPDATE").
         */
        void File(const std::string& filename, const std::string& option = "RECREATE");

    private:
        KinematicsProcessor& _proc;
        ConfigReaction& _rad; 
        PhysicsSelection* _sel = nullptr; 
        std::string _maskCol;
        bool _initialized = false;
        
        // Style Guide Fix: Internal collections use std::vector, not RVec
        std::vector<SplitDef> _splits;
        std::vector<HistoDef> _defs; 
        std::vector<HistoDef2D> _defs2D; 
        
        std::vector<ROOT::RDF::RResultPtr<THnSparseD>> _results;

        /** @brief Internal helper to book a 1D histogram. */
        void BookInternal(const HistoDef& def);

        /** @brief Internal helper to book a 2D histogram. */
        void BookInternal2D(const HistoDef2D& def);

        /** @brief Recursive helper to project N-D histogram to 1D slices. */
        void UnfoldAndWrite(THnSparseD* hn, TDirectory* dir, std::string current_suffix = "", int axis_depth = 1);

        /** @brief Safely checks if a column exists and throws a descriptive error. */
        void CheckColumn(const std::string& baseName, const std::string& fullName);
    };

    // =================================================================================
    // IMPLEMENTATION: Histogrammer
    // =================================================================================

    inline Histogrammer::Histogrammer(KinematicsProcessor& proc, PhysicsSelection* sel) 
        : _proc(proc), _rad(*proc.Reaction()), _maskCol("") 
    {
        if (sel) _sel = sel;
    }

    inline void Histogrammer::AddSplit(const std::string& name, const std::string& baseCol, int nbins, double min, double max) {
        _splits.push_back({name, baseCol, nbins, min, max});
    }

    inline void Histogrammer::Create(const std::string& name, const std::string& title, int nbins, double min, double max, const std::string& varBaseName) {
        _defs.push_back({name, title, nbins, min, max, varBaseName});
    }
    
    inline void Histogrammer::Create2D(const std::string& name, const std::string& title, int nbinsX, double xmin, double xmax, int nbinsY, double ymin, double ymax, const std::string& XvarBaseName, const std::string& YvarBaseName) {
        _defs2D.push_back({name, title, nbinsX, xmin, xmax, nbinsY, ymin, ymax, XvarBaseName, YvarBaseName});
    }

    inline void Histogrammer::CheckColumn(const std::string& baseName, const std::string& fullName) {
        // Style Guide: Safe Getter Pattern
        if (!_rad.ColumnExists(fullName)) {
            std::string err = "\n\n[RAD Histogrammer ERROR] Missing Column: '" + fullName + "'\n";
            err += " -> Stream Prefix : '" + _proc.GetPrefix() + "'\n";
            err += " -> Requested Base: '" + baseName + "'\n";
            
            if (baseName.find("rec_") == 0 || baseName.find("tru_") == 0 || baseName.find("mc_") == 0) {
                err += "\n[!] CAUSE: You included the stream prefix in your histogram recipe.";
                err += "\n    FIX: Use the base name only. The framework handles prefixes automatically.\n\n";
            } 
            else if (_proc.GetPrefix().find("tru") != std::string::npos) {
                err += "\n[!] CAUSE: You are trying to plot a detector variable while processing a Truth stream.";
                err += "\n    FIX: Detector variables only exist in the Rec stream. Separate your detector plots into a recipe specifically for Rec().\n\n";
            } 
            else {
                err += "\n[!] CAUSE: Variable does not exist. Check spelling or ensure the calculation is registered in your topology_recipe.\n\n";
            }
            throw std::runtime_error(err);
        }
    }
    
    inline void Histogrammer::Init() {
        if (_initialized) return;
        if (_sel) _maskCol = _sel->GetMaskColumn();
        for (const auto& def : _defs) BookInternal(def);
        for (const auto& def2 : _defs2D) BookInternal2D(def2);
        _initialized = true;
    }
    
    inline void Histogrammer::BookInternal(const HistoDef& def) {
        std::string fullVarName = _proc.FullName(def.varBaseName);
        CheckColumn(def.varBaseName, fullVarName);
        
        // Style Guide: Column names as std::vector for RDF Interface
        std::vector<std::string> cols;
        cols.push_back(fullVarName);

        int n_dims = 1 + _splits.size();
        std::vector<int> bins_vec(n_dims);
        std::vector<double> xmin(n_dims), xmax(n_dims);
        
        bins_vec[0] = def.nbins; xmin[0] = def.min; xmax[0] = def.max;
        
        for(size_t i=0; i<_splits.size(); ++i) {
            std::string fullSplitCol = _proc.FullName(_splits[i].baseCol);
            CheckColumn(_splits[i].baseCol, fullSplitCol);
            
            cols.push_back(fullSplitCol);
            bins_vec[i+1] = _splits[i].nbins;
            xmin[i+1] = _splits[i].min;
            xmax[i+1] = _splits[i].max;
        }

        std::string hNameFull = _proc.GetPrefix() + def.name + _proc.GetSuffix();
        auto hist = std::make_shared<THnSparseD>(hNameFull.c_str(), def.title.c_str(), n_dims, bins_vec.data(), xmin.data(), xmax.data());
        
        hist->GetAxis(0)->SetName(def.name.c_str());
        for(size_t i=0; i<_splits.size(); ++i) hist->GetAxis(i+1)->SetName(_splits[i].name.c_str());

        if(!_maskCol.empty()) {
            cols.push_back(_maskCol); 
            THnCombi<true> action(hist);
            _results.push_back(_rad.CurrFrame().Book(std::move(action), cols));
        } else {
            THnCombi<false> action(hist);
            _results.push_back(_rad.CurrFrame().Book(std::move(action), cols));
        }
    }
    
    inline void Histogrammer::BookInternal2D(const HistoDef2D& def) {
        std::string xName = _proc.FullName(def.XvarBaseName);
        std::string yName = _proc.FullName(def.YvarBaseName);

        CheckColumn(def.XvarBaseName, xName);
        CheckColumn(def.YvarBaseName, yName);

        std::vector<std::string> cols = {xName, yName};

        int n_dims = 2 + _splits.size();
        std::vector<int> bins_vec(n_dims);
        std::vector<double> xmin(n_dims), xmax(n_dims);

        bins_vec[0] = def.nbinsX; xmin[0] = def.xmin; xmax[0] = def.xmax;
        bins_vec[1] = def.nbinsY; xmin[1] = def.ymin; xmax[1] = def.ymax;

        for (size_t i = 0; i < _splits.size(); ++i) {
            auto& s = _splits[i];
            std::string col = _proc.FullName(s.baseCol);
            
            CheckColumn(s.baseCol, col);
            
            cols.push_back(col);
            bins_vec[i + 2] = s.nbins;
            xmin[i + 2] = s.min;
            xmax[i + 2] = s.max;
        }

        std::string fullName = _proc.GetPrefix() + def.name + _proc.GetSuffix();
        auto hist = std::make_shared<THnSparseD>(fullName.c_str(), def.title.c_str(), n_dims, bins_vec.data(), xmin.data(), xmax.data());

        hist->GetAxis(0)->SetName(def.XvarBaseName.c_str());
        hist->GetAxis(1)->SetName(def.YvarBaseName.c_str());
        for (size_t i = 0; i < _splits.size(); ++i) hist->GetAxis(i + 2)->SetName(_splits[i].name.c_str());

        if (!_maskCol.empty()) {
            cols.push_back(_maskCol);
            THnCombi<true> action(hist);
            _results.push_back(_rad.CurrFrame().Book(std::move(action), cols));
        } else {
            THnCombi<false> action(hist);
            _results.push_back(_rad.CurrFrame().Book(std::move(action), cols));
        }
    }
    
    inline void Histogrammer::File(const std::string& filename, const std::string& option) {
        auto file = std::unique_ptr<TFile>(TFile::Open(filename.c_str(), option.c_str())); 
        if (!file || file->IsZombie()) {
            std::cerr << "[RAD Histogrammer ERROR] Could not open file: " << filename << std::endl;
            return;
        }
        for(auto& resultPtr : _results) {
            auto hist_raw = resultPtr->Clone();
            auto hn = dynamic_cast<THnSparseD*>(hist_raw);
            if (hn) {
                TDirectory* subdir = file->mkdir(hn->GetName());
                if (!subdir) subdir = file->GetDirectory(hn->GetName());
                subdir->cd();
                UnfoldAndWrite(hn, subdir);
            }
            delete hn; 
        }
        file->Close();
    }

    inline void Histogrammer::UnfoldAndWrite(THnSparseD* hn, TDirectory* dir, std::string current_suffix, int axis_depth) {
        int ndim = hn->GetNdimensions();
        int nsplits = _splits.size();
        bool is2D = (ndim == nsplits + 2);

        if (axis_depth > nsplits) {
            dir->cd();
            TH1* proj = nullptr;
            if (is2D) proj = hn->Projection(0, 1, "E");
            else      proj = hn->Projection(0, "E");

            std::string full_name = std::string(hn->GetName()) + current_suffix;
            proj->SetName(full_name.c_str());
            proj->SetTitle((std::string(hn->GetTitle()) + current_suffix).c_str());

            proj->Write();
            delete proj;
            return;
        }

        auto axis = hn->GetAxis(axis_depth);
        int nbins = axis->GetNbins();
        int saved_min = axis->GetFirst();
        int saved_max = axis->GetLast();

        for (int i = 1; i <= nbins; ++i) {
            axis->SetRange(i, i);
            std::string suffix = current_suffix + "_" + axis->GetName() + std::to_string(i);
            UnfoldAndWrite(hn, dir, suffix, axis_depth + 1);
        }
        axis->SetRange(saved_min, saved_max);
    }

} // namespace histo
} // namespace rad
