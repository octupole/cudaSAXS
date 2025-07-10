#include "AdvancedOptionsDialog.h"
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QLabel>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QLineEdit>
#include <QComboBox>
#include <QCheckBox>
#include <QDialogButtonBox>
#include <QStackedWidget>

AdvancedOptionsDialog::AdvancedOptionsDialog(QWidget *parent)
    : QDialog(parent)
{
    setWindowTitle("Advanced Options");
    setModal(true);
    setupUI();
    connectSignals();
}

void AdvancedOptionsDialog::setupUI()
{
    auto mainLayout = new QVBoxLayout(this);

    // B-spline order
    auto bsplineGroup = new QGroupBox("B-Spline Settings");
    auto bsplineLayout = new QHBoxLayout(bsplineGroup);
    bsplineLayout->addWidget(new QLabel("B-Spline Order:"));
    m_orderSpin = new QSpinBox;
    m_orderSpin->setRange(1, 10);
    m_orderSpin->setValue(4);
    bsplineLayout->addWidget(m_orderSpin);
    bsplineLayout->addStretch();
    mainLayout->addWidget(bsplineGroup);

    // Scaled grid settings
    auto scaledGroup = new QGroupBox("Scaled Grid Settings");
    auto scaledLayout = new QVBoxLayout(scaledGroup);

    m_scaledGridCheck = new QCheckBox("Use custom scaled grid");
    scaledLayout->addWidget(m_scaledGridCheck);

    m_scaledGridWidget = new QWidget;
    auto scaledGridLayout = new QVBoxLayout(m_scaledGridWidget);

    // Scale factor
    auto scaleLayout = new QHBoxLayout;
    scaleLayout->addWidget(new QLabel("Scale Factor (σ):"));
    m_scaleFactorSpin = new QDoubleSpinBox;
    m_scaleFactorSpin->setRange(0.1, 10.0);
    m_scaleFactorSpin->setSingleStep(0.1);
    m_scaleFactorSpin->setValue(2.5);
    m_scaleFactorSpin->setDecimals(1);
    scaleLayout->addWidget(m_scaleFactorSpin);
    scaleLayout->addStretch();
    scaledGridLayout->addLayout(scaleLayout);

    // Scaled grid type
    auto scaledTypeLayout = new QHBoxLayout;
    scaledTypeLayout->addWidget(new QLabel("Scaled Grid Type:"));
    m_scaledGridTypeCombo = new QComboBox;
    m_scaledGridTypeCombo->addItems({"Cubic (single value)", "Custom (x, y, z)"});
    scaledTypeLayout->addWidget(m_scaledGridTypeCombo);
    scaledTypeLayout->addStretch();
    scaledGridLayout->addLayout(scaledTypeLayout);

    // Stacked widget for scaled grid input
    auto scaledStack = new QStackedWidget;

    // Single scaled grid value
    m_singleScaledWidget = new QWidget;
    auto singleScaledLayout = new QHBoxLayout(m_singleScaledWidget);
    singleScaledLayout->addWidget(new QLabel("Scaled Grid Size:"));
    m_scaledGridSingleSpin = new QSpinBox;
    m_scaledGridSingleSpin->setRange(1, 2000);
    m_scaledGridSingleSpin->setValue(160);
    singleScaledLayout->addWidget(m_scaledGridSingleSpin);
    singleScaledLayout->addStretch();

    // Triple scaled grid values
    m_tripleScaledWidget = new QWidget;
    auto tripleScaledLayout = new QHBoxLayout(m_tripleScaledWidget);
    tripleScaledLayout->addWidget(new QLabel("X:"));
    m_scaledGridXSpin = new QSpinBox;
    m_scaledGridXSpin->setRange(1, 2000);
    m_scaledGridXSpin->setValue(160);
    tripleScaledLayout->addWidget(m_scaledGridXSpin);

    tripleScaledLayout->addWidget(new QLabel("Y:"));
    m_scaledGridYSpin = new QSpinBox;
    m_scaledGridYSpin->setRange(1, 2000);
    m_scaledGridYSpin->setValue(160);
    tripleScaledLayout->addWidget(m_scaledGridYSpin);

    tripleScaledLayout->addWidget(new QLabel("Z:"));
    m_scaledGridZSpin = new QSpinBox;
    m_scaledGridZSpin->setRange(1, 2000);
    m_scaledGridZSpin->setValue(160);
    tripleScaledLayout->addWidget(m_scaledGridZSpin);
    tripleScaledLayout->addStretch();

    scaledStack->addWidget(m_singleScaledWidget);
    scaledStack->addWidget(m_tripleScaledWidget);
    scaledGridLayout->addWidget(scaledStack);

    connect(m_scaledGridTypeCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
            scaledStack, &QStackedWidget::setCurrentIndex);

    m_scaledGridWidget->setEnabled(false);
    scaledLayout->addWidget(m_scaledGridWidget);

    mainLayout->addWidget(scaledGroup);

    // Water model settings
    auto waterGroup = new QGroupBox("Water Model / Padding");
    auto waterLayout = new QVBoxLayout(waterGroup);

    m_waterModelCheck = new QCheckBox("Use specific water model");
    waterLayout->addWidget(m_waterModelCheck);

    m_waterModelWidget = new QWidget;
    auto waterModelLayout = new QGridLayout(m_waterModelWidget);

    waterModelLayout->addWidget(new QLabel("Water Model:"), 0, 0);
    m_waterModelEdit = new QLineEdit;
    m_waterModelEdit->setPlaceholderText("e.g., SPCE, TIP3P");
    waterModelLayout->addWidget(m_waterModelEdit, 0, 1);

    waterModelLayout->addWidget(new QLabel("Sodium Atoms:"), 1, 0);
    m_sodiumSpin = new QSpinBox;
    m_sodiumSpin->setRange(0, 10000);
    waterModelLayout->addWidget(m_sodiumSpin, 1, 1);

    waterModelLayout->addWidget(new QLabel("Chlorine Atoms:"), 2, 0);
    m_chlorineSpin = new QSpinBox;
    m_chlorineSpin->setRange(0, 10000);
    waterModelLayout->addWidget(m_chlorineSpin, 2, 1);

    m_waterModelWidget->setEnabled(false);
    waterLayout->addWidget(m_waterModelWidget);

    mainLayout->addWidget(waterGroup);

    // Dialog buttons
    auto buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    connect(buttonBox, &QDialogButtonBox::accepted, this, &QDialog::accept);
    connect(buttonBox, &QDialogButtonBox::rejected, this, &QDialog::reject);
    mainLayout->addWidget(buttonBox);

    mainLayout->addStretch();
}

void AdvancedOptionsDialog::connectSignals()
{
    connect(m_scaledGridCheck, &QCheckBox::toggled,
            this, &AdvancedOptionsDialog::onScaledGridEnabledChanged);
    connect(m_scaledGridTypeCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
            this, &AdvancedOptionsDialog::onScaledGridTypeChanged);
    connect(m_waterModelCheck, &QCheckBox::toggled,
            this, &AdvancedOptionsDialog::onWaterModelEnabledChanged);
}

void AdvancedOptionsDialog::onScaledGridEnabledChanged(bool enabled)
{
    m_scaledGridWidget->setEnabled(enabled);
}

void AdvancedOptionsDialog::onScaledGridTypeChanged(int index)
{
    // Nothing special needed here as the stacked widget handles it
}

void AdvancedOptionsDialog::onWaterModelEnabledChanged(bool enabled)
{
    m_waterModelWidget->setEnabled(enabled);
}

void AdvancedOptionsDialog::setOptionsData(const OptionsData &data)
{
    m_orderSpin->setValue(data.bsplineOrder);
    m_scaleFactorSpin->setValue(data.scaleFactor);

    // Scaled grid
    bool hasScaledGrid = !data.scaledGrid.isEmpty();
    m_scaledGridCheck->setChecked(hasScaledGrid);

    if (hasScaledGrid)
    {
        if (data.scaledGrid.size() == 1)
        {
            m_scaledGridTypeCombo->setCurrentIndex(0);
            m_scaledGridSingleSpin->setValue(data.scaledGrid[0]);
        }
        else if (data.scaledGrid.size() == 3)
        {
            m_scaledGridTypeCombo->setCurrentIndex(1);
            m_scaledGridXSpin->setValue(data.scaledGrid[0]);
            m_scaledGridYSpin->setValue(data.scaledGrid[1]);
            m_scaledGridZSpin->setValue(data.scaledGrid[2]);
        }
    }

    // Water model
    bool hasWaterModel = !data.waterModel.isEmpty();
    m_waterModelCheck->setChecked(hasWaterModel);

    if (hasWaterModel)
    {
        m_waterModelEdit->setText(data.waterModel);
        m_sodiumSpin->setValue(data.sodiumAtoms);
        m_chlorineSpin->setValue(data.chlorineAtoms);
    }
}

OptionsData AdvancedOptionsDialog::getOptionsData() const
{
    OptionsData data;

    data.bsplineOrder = m_orderSpin->value();
    data.scaleFactor = m_scaleFactorSpin->value();

    // Scaled grid
    if (m_scaledGridCheck->isChecked())
    {
        data.scaledGrid.clear();
        if (m_scaledGridTypeCombo->currentIndex() == 0)
        {
            data.scaledGrid.append(m_scaledGridSingleSpin->value());
        }
        else
        {
            data.scaledGrid.append(m_scaledGridXSpin->value());
            data.scaledGrid.append(m_scaledGridYSpin->value());
            data.scaledGrid.append(m_scaledGridZSpin->value());
        }
    }

    // Water model
    if (m_waterModelCheck->isChecked())
    {
        data.waterModel = m_waterModelEdit->text();
        data.sodiumAtoms = m_sodiumSpin->value();
        data.chlorineAtoms = m_chlorineSpin->value();
    }

    return data;
}