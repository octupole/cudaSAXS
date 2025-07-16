#include "InputForm.h"
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QComboBox>
#include <QPushButton>
#include <QFileDialog>
#include <QStackedWidget>

InputForm::InputForm(QWidget *parent) : QWidget(parent)
{
    setupUI();
    connectSignals();
    resetToDefaults();
    validateForm();
}

void InputForm::setupUI()
{
    auto mainLayout = new QVBoxLayout(this);

    // File selection group
    auto fileGroup = new QGroupBox("Input Files");
    auto fileLayout = new QGridLayout(fileGroup);

    fileLayout->addWidget(new QLabel("Topology (.tpr):"), 0, 0);
    m_tprFileEdit = new QLineEdit;
    m_tprBrowseBtn = new QPushButton("Browse...");
    fileLayout->addWidget(m_tprFileEdit, 0, 1);
    fileLayout->addWidget(m_tprBrowseBtn, 0, 2);

    fileLayout->addWidget(new QLabel("Trajectory (.xtc):"), 1, 0);
    m_xtcFileEdit = new QLineEdit;
    m_xtcBrowseBtn = new QPushButton("Browse...");
    fileLayout->addWidget(m_xtcFileEdit, 1, 1);
    fileLayout->addWidget(m_xtcBrowseBtn, 1, 2);

    mainLayout->addWidget(fileGroup);

    // Grid settings group
    auto gridGroup = new QGroupBox("Grid Settings");
    auto gridLayout = new QVBoxLayout(gridGroup);

    auto gridTypeLayout = new QHBoxLayout;
    gridTypeLayout->addWidget(new QLabel("Grid Type:"));
    m_gridTypeCombo = new QComboBox;
    m_gridTypeCombo->addItems({"Cubic (single value)", "Custom (x, y, z)"});
    gridTypeLayout->addWidget(m_gridTypeCombo);
    gridTypeLayout->addStretch();
    gridLayout->addLayout(gridTypeLayout);

    // Stacked widget for grid input
    auto gridStack = new QStackedWidget;

    // Single grid value
    m_singleGridWidget = new QWidget;
    auto singleLayout = new QHBoxLayout(m_singleGridWidget);
    singleLayout->addWidget(new QLabel("Grid Size:"));
    m_gridSingleSpin = new QSpinBox;
    m_gridSingleSpin->setRange(1, 1000);
    m_gridSingleSpin->setValue(64);
    singleLayout->addWidget(m_gridSingleSpin);
    singleLayout->addStretch();

    // Triple grid values
    m_tripleGridWidget = new QWidget;
    auto tripleLayout = new QHBoxLayout(m_tripleGridWidget);
    tripleLayout->addWidget(new QLabel("X:"));
    m_gridXSpin = new QSpinBox;
    m_gridXSpin->setRange(1, 1000);
    m_gridXSpin->setValue(64);
    tripleLayout->addWidget(m_gridXSpin);

    tripleLayout->addWidget(new QLabel("Y:"));
    m_gridYSpin = new QSpinBox;
    m_gridYSpin->setRange(1, 1000);
    m_gridYSpin->setValue(64);
    tripleLayout->addWidget(m_gridYSpin);

    tripleLayout->addWidget(new QLabel("Z:"));
    m_gridZSpin = new QSpinBox;
    m_gridZSpin->setRange(1, 1000);
    m_gridZSpin->setValue(64);
    tripleLayout->addWidget(m_gridZSpin);
    tripleLayout->addStretch();

    gridStack->addWidget(m_singleGridWidget);
    gridStack->addWidget(m_tripleGridWidget);
    gridLayout->addWidget(gridStack);

    connect(m_gridTypeCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
            gridStack, &QStackedWidget::setCurrentIndex);

    mainLayout->addWidget(gridGroup);

    // Frame settings group
    auto frameGroup = new QGroupBox("Frame Settings");
    auto frameLayout = new QGridLayout(frameGroup);

    frameLayout->addWidget(new QLabel("Start Frame:"), 0, 0);
    m_startFrameSpin = new QSpinBox;
    m_startFrameSpin->setRange(0, 999999);
    frameLayout->addWidget(m_startFrameSpin, 0, 1);

    frameLayout->addWidget(new QLabel("End Frame:"), 0, 2);
    m_endFrameSpin = new QSpinBox;
    m_endFrameSpin->setRange(0, 999999);
    m_endFrameSpin->setValue(100);
    frameLayout->addWidget(m_endFrameSpin, 0, 3);

    frameLayout->addWidget(new QLabel("Frame Interval:"), 1, 0);
    m_frameIntervalSpin = new QSpinBox;
    m_frameIntervalSpin->setRange(1, 1000);
    m_frameIntervalSpin->setValue(1);
    frameLayout->addWidget(m_frameIntervalSpin, 1, 1);

    mainLayout->addWidget(frameGroup);

    // Basic options group
    auto optionsGroup = new QGroupBox("Basic Options");
    auto optionsLayout = new QGridLayout(optionsGroup);

    optionsLayout->addWidget(new QLabel("Output File:"), 0, 0);
    m_outputFileEdit = new QLineEdit("saxs.dat");
    m_outputBrowseBtn = new QPushButton("Browse...");
    optionsLayout->addWidget(m_outputFileEdit, 0, 1);
    optionsLayout->addWidget(m_outputBrowseBtn, 0, 2);

    optionsLayout->addWidget(new QLabel("Bin Size (Dq):"), 1, 0);
    m_binSizeSpin = new QDoubleSpinBox;
    m_binSizeSpin->setRange(0.001, 1.0);
    m_binSizeSpin->setSingleStep(0.01);
    m_binSizeSpin->setValue(0.05);
    m_binSizeSpin->setDecimals(3);
    optionsLayout->addWidget(m_binSizeSpin, 1, 1);

    optionsLayout->addWidget(new QLabel("Q Cutoff:"), 2, 0);
    m_qCutoffSpin = new QDoubleSpinBox;
    m_qCutoffSpin->setRange(0.1, 10.0);
    m_qCutoffSpin->setSingleStep(0.1);
    m_qCutoffSpin->setValue(4.0);
    m_qCutoffSpin->setDecimals(1);
    optionsLayout->addWidget(m_qCutoffSpin, 2, 1);

    optionsLayout->addWidget(new QLabel("Simulation Type:"), 3, 0);
    m_simulationTypeCombo = new QComboBox;
    m_simulationTypeCombo->addItems({"npt", "nvt"});
    optionsLayout->addWidget(m_simulationTypeCombo, 3, 1);

    mainLayout->addWidget(optionsGroup);

    // Advanced options button
    m_advancedBtn = new QPushButton("Advanced Options...");
    mainLayout->addWidget(m_advancedBtn);

    // Status label
    m_statusLabel = new QLabel;
    m_statusLabel->setWordWrap(true);
    mainLayout->addWidget(m_statusLabel);

    // Run button
    m_runBtn = new QPushButton("Run SAXS Calculation");
    m_runBtn->setDefault(true);
    mainLayout->addWidget(m_runBtn);

    mainLayout->addStretch();
}

void InputForm::connectSignals()
{
    connect(m_tprBrowseBtn, &QPushButton::clicked, this, &InputForm::onBrowseTPR);
    connect(m_xtcBrowseBtn, &QPushButton::clicked, this, &InputForm::onBrowseXTC);
    connect(m_outputBrowseBtn, &QPushButton::clicked, this, &InputForm::onBrowseOutput);
    connect(m_gridTypeCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
            this, &InputForm::onGridTypeChanged);

    connect(m_advancedBtn, &QPushButton::clicked, this, &InputForm::advancedOptionsRequested);
    connect(m_runBtn, &QPushButton::clicked, this, &InputForm::runRequested);

    // Connect all inputs to validation
    connect(m_tprFileEdit, &QLineEdit::textChanged, this, &InputForm::validateForm);
    connect(m_xtcFileEdit, &QLineEdit::textChanged, this, &InputForm::validateForm);
    connect(m_gridSingleSpin, QOverload<int>::of(&QSpinBox::valueChanged), this, &InputForm::validateForm);
    connect(m_gridXSpin, QOverload<int>::of(&QSpinBox::valueChanged), this, &InputForm::validateForm);
    connect(m_gridYSpin, QOverload<int>::of(&QSpinBox::valueChanged), this, &InputForm::validateForm);
    connect(m_gridZSpin, QOverload<int>::of(&QSpinBox::valueChanged), this, &InputForm::validateForm);
    connect(m_startFrameSpin, QOverload<int>::of(&QSpinBox::valueChanged), this, &InputForm::validateForm);
    connect(m_endFrameSpin, QOverload<int>::of(&QSpinBox::valueChanged), this, &InputForm::validateForm);
}

void InputForm::onBrowseTPR()
{
    QString filename = QFileDialog::getOpenFileName(this,
                                                    "Select Topology File", QString(), "TPR Files (*.tpr);;All Files (*)");
    if (!filename.isEmpty())
    {
        m_tprFileEdit->setText(filename);
    }
}

void InputForm::onBrowseXTC()
{
    QString filename = QFileDialog::getOpenFileName(this,
                                                    "Select Trajectory File", QString(), "XTC Files (*.xtc);;All Files (*)");
    if (!filename.isEmpty())
    {
        m_xtcFileEdit->setText(filename);
    }
}

void InputForm::onBrowseOutput()
{
    QString filename = QFileDialog::getSaveFileName(this,
                                                    "Select Output File", QString(), "Data Files (*.dat);;All Files (*)");
    if (!filename.isEmpty())
    {
        m_outputFileEdit->setText(filename);
    }
}

void InputForm::onGridTypeChanged(int index)
{
    validateForm();
}

void InputForm::validateForm()
{
    OptionsData data = getOptionsData();
    bool isValid = data.isValid();
    QString error = data.validationError();

    m_runBtn->setEnabled(isValid);

    if (isValid)
    {
        m_statusLabel->setText("Ready to run");
        m_statusLabel->setStyleSheet("QLabel { color: green; }");
    }
    else
    {
        m_statusLabel->setText(error);
        m_statusLabel->setStyleSheet("QLabel { color: red; }");
    }

    emit formValidityChanged(isValid);
}

OptionsData InputForm::getOptionsData() const
{
    OptionsData data = m_advancedOptions; // Start with advanced options

    // Override with form values
    data.tprFile = m_tprFileEdit->text();
    data.xtcFile = m_xtcFileEdit->text();

    // Grid size
    data.gridSize.clear();
    if (m_gridTypeCombo->currentIndex() == 0)
    {
        data.gridSize.append(m_gridSingleSpin->value());
    }
    else
    {
        data.gridSize.append(m_gridXSpin->value());
        data.gridSize.append(m_gridYSpin->value());
        data.gridSize.append(m_gridZSpin->value());
    }

    data.startFrame = m_startFrameSpin->value();
    data.endFrame = m_endFrameSpin->value();
    data.frameInterval = m_frameIntervalSpin->value();

    data.outputFile = m_outputFileEdit->text();
    data.binSize = m_binSizeSpin->value();
    data.qCutoff = m_qCutoffSpin->value();
    data.simulationType = m_simulationTypeCombo->currentText();

    return data;
}

void InputForm::setOptionsData(const OptionsData &data)
{
    m_tprFileEdit->setText(data.tprFile);
    m_xtcFileEdit->setText(data.xtcFile);

    // Grid size
    if (data.gridSize.size() == 1)
    {
        m_gridTypeCombo->setCurrentIndex(0);
        m_gridSingleSpin->setValue(data.gridSize[0]);
    }
    else if (data.gridSize.size() == 3)
    {
        m_gridTypeCombo->setCurrentIndex(1);
        m_gridXSpin->setValue(data.gridSize[0]);
        m_gridYSpin->setValue(data.gridSize[1]);
        m_gridZSpin->setValue(data.gridSize[2]);
    }

    m_startFrameSpin->setValue(data.startFrame);
    m_endFrameSpin->setValue(data.endFrame);
    m_frameIntervalSpin->setValue(data.frameInterval);

    m_outputFileEdit->setText(data.outputFile);
    m_binSizeSpin->setValue(data.binSize);
    m_qCutoffSpin->setValue(data.qCutoff);

    int simIndex = m_simulationTypeCombo->findText(data.simulationType);
    if (simIndex >= 0)
    {
        m_simulationTypeCombo->setCurrentIndex(simIndex);
    }

    m_advancedOptions = data;
    validateForm();
}

void InputForm::resetToDefaults()
{
    OptionsData defaults;
    defaults.gridSize.append(64);
    defaults.startFrame = 0;
    defaults.endFrame = 100;
    setOptionsData(defaults);
}